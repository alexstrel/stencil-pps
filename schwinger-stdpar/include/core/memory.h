#pragma once
#include <common.h>
#include <memory>
#include <memory_resource>

constexpr std::size_t alignment_req = std::alignment_of_v<long double>;

constexpr int max_n_buffers = 128;

class UpstreamMemoryResource : public std::pmr::memory_resource {
  protected:
    //
    std::size_t nBytes = 0;
    
    void* do_allocate(std::size_t nbytes, std::size_t alignment = alignment_req) override {//alignof(std::max_align)
      nBytes = ( (nbytes + alignment - 1) / alignment ) *  alignment;
      //
      return std::aligned_alloc(alignment, nBytes);
    }
    
    void do_deallocate(void* ptr, std::size_t nbytes, std::size_t alignment = alignment_req) override {
      //
      return std::free(ptr);
    }
    
    bool do_is_equal(const std::pmr::memory_resource &other) const noexcept override { return (this == &other); } 
    
  public:
  
    UpstreamMemoryResource()                                = default;
    UpstreamMemoryResource(const UpstreamMemoryResource &)  = default;
    UpstreamMemoryResource(UpstreamMemoryResource &&)       = default;     
};

static auto upstream = UpstreamMemoryResource{};

class PMRBuffer{
  public:
    //
    std::size_t                          pmr_base_bytes;
    std::size_t                          pmr_bytes;
    std::shared_ptr<std::byte[]>         pmr_ptr;//needed only for pmr-used spinor    
    
    std::shared_ptr<std::pmr::monotonic_buffer_resource>  pmr_pool;    

    PMRState pmr_state = PMRState::InvalidState; 
    //
    bool is_exclusive;

    PMRBuffer(const std::size_t nbytes, const PMRState state = PMRState::Vacant) : pmr_base_bytes(nbytes),
                                                                                   pmr_bytes(((nbytes + alignment_req - 1) / alignment_req ) *  alignment_req), 
                                                                                   pmr_ptr(std::make_shared<std::byte[]>(pmr_bytes)),
                                                                                   pmr_pool(std::make_shared<std::pmr::monotonic_buffer_resource>(pmr_ptr.get(), pmr_bytes)),                          
                                                                                   pmr_state(state),
                                                                                   is_exclusive( state != PMRState::Locked and state != PMRState::Reserved ? false : true) { } 

    PMRBuffer()                             = default;
    PMRBuffer(const PMRBuffer &)            = default;
    PMRBuffer(PMRBuffer &&)                 = default;

    virtual ~PMRBuffer()                    = default;

    auto Data()  const { return pmr_ptr;       }
    auto Get()   const { return pmr_ptr.get(); }
    
    auto BaseBytes()  const { return pmr_base_bytes; }        
    
    auto Bytes() const { return pmr_bytes; }    
    auto State() const { return pmr_state; }
    //
    auto Pool()  const { return pmr_pool;  } 
    
    bool IsExclusive const { return is_exclusive; }           
    
    void SetState(PMRState state) { 
      if ((pmr_state == PMRState::Vacant) and ( state == PMRState::Locked or state == PMRState::Reserved)) is_exclusive = true;
      pmr_state = state; 
    }
    
    void ResetPool() const {  pmr_pool->reset(); } 
    
    bool IsReserved(const std::size_t nbytes) const {  
      if (pmr_ptr != nullptr ) {
        if (pmr_state != PMRState::Reserved) {
          return false;
        }
        //
        return (nbytes <= base_bytes);
      } 
      //
      return false;
    }    
    
    void Release() {
      if (is_exclusive) {
        pmr_state = PMRState::Vacant;
        //
      } else {
        if (pmr_state == PMRState::Locked) pmr_state = PMRState::NonVacant;
        else                               pmr_state = PMRState::Vacant;
      }
      pmr_pool->reset();
      //      
      is_exclusive = false;
    }
               
    void Destroy()   const {  
      pmr_pool.reset(nullptr);
      pmr_base_bytes = 0;
      pmr_bytes      = 0;      
      pmr_ptr.reset(nullptr);
      //
      pmr_state = PMRState::InvalidState;
    } 
    
    bool IsEqual(const PMRBuffer &other) const { return (this->pmr_ptr.get() == other.Get()); }               
};

static PMRBuffer default_pmr_buffer = PMRBuffer();

/**
    Derived from the QUDA library, original developer : Kate Clark (NVIDIA)
*/

using PMRBufferRef = std::reference_wrapper<PMRBuffer>;

namespace pmr_pool
{
  /**
      Contains all instanses of pmr buffer (incl. unallocated)
  */
  static std::vector<std::shared_ptr<PMRBuffer>> pmrBuffers;
  
  /** Cache of pmr-memory allocations.  We cache pmr
      memory allocations so that fields can reuse them.*/
  static std::multimap<std::size_t, std::shared_ptr<PMRBuffer>> pmrCachedBuffers;

  /** Sizes of active pmr-memory allocations.  For convenience,
      we keep track of the sizes of active allocations (i.e., those not
      in the cache). */
  static std::map<std::byte*, std::shared_ptr<PMRBuffer>> pmrAllocatedBuffers;
  
  static bool pmr_memory_pool = false;   

  template<bool is_exclusive = true>
  decltype(auto) pmr_malloc(const std::size_t nbytes, const bool is_reserved = false) {
  
    PRMState final_buffer_state = is_reserved ? PMRState::Reserved : (is_exclusive ? PMRState::Locked : PMRState::NonVacant);    
     
    if( pmr_memory_pool == false ) { 
      pmrBuffers.reserve(max_n_buffers);
      //
      pmr_memory_pool = true;
    }
        
    if ( not pmrCachedBuffers.empty() ) {//if buffer is not empty... 

      auto buffer_range = pmrCachedBuffers.equal_range(nbytes);      
          
      if ( buffer_range.first != pmrCachedBuffers.end() ) {//...and threre is one (not locked or reserved) with sufficient size
        //
        auto subrange_view = std::range::subrange(buffer_range.first, buffer_range.second);
	  
        for(auto& [bytes, pmr_buffer] : subrange_view) {
	  if ( ((pmr_buffer->State() == PMRState::Vacant) or (pmr_buffer->State() == PMRState::NonVacant and not is_exclusive)) and (pmr_buffer->State() != PMRState::Reserved) ) {
	    //
            if ( pmr_buffer->State() == PMRState::NonVacant ) { //this is a shared buffer
              pmr_buffer->ResetPool(); 
              // we locked buffer if it's not reserved:
              if (not is_reserved) final_buffer_state = PMRState::Locked; 
            }
	    //
	    pmr_buffer->SetState( final_buffer_state );
	    //	    
	    return pmr_buffer;
	  } 
        } 
      }
    }
    //
    //auto pmr_buffer = PMRBuffer{nbytes, buffer_state};
    pmrBuffers.emplace_back(std::make_shared<PMRBuffer>(nbytes, final_buffer_state));
    //
    auto pmr_buffer = pmrBuffers.back();
    //
    pmrCachedBuffers.insert(std::make_pair(nbytes, pmr_buffer));	
    //
    pmrAllocatedBuffers.insert(std::make_pair(pmr_buffer->Get(), pmr_buffer));
        
    return pmr_buffer;
  }    

  template<ArithmeticTp T>
  void pmr_free(T* ptr) {
    //
    auto raw_ptr = reinterpret_cast<std::byte*>(ptr);
    
    if ( pmrAllocatedBuffers.empty() or (not pmrAllocatedBuffers.contains(raw_ptr)) ) { std::cerr << "Attempt to free invalid pmr buffer" << std::endl; }    
    
    auto it = pmrAllocatedBuffers.find(raw_ptr);
    
    if (it != pmrAllocatedBuffers.end()) {
    
      auto& pmr_buffer   = it->second; 
      //
      if (pmr_buffer->State() != PMRState::Vacant) { std::cerr << "Attempt to free pmr buffer that is still in use." << std::endl; }           
      //
      std::size_t nbytes = pmr_buffer->BaseBytes();
      //
      auto range = pmrCachedBuffers.equal_range(nbytes);      
      //
      for (auto& range_it = range.first; range_it != range.second; ++range_it) {
        auto &cached_pmr_buffer = range_it->second;
        if ( cached_pmr_buffer->IsEquel(pmr_buffer) ) {
          //
          pmrAllocatedBuffers.erase(range_it);                    
          //
          break;
        }
      }     
      //
      pmr_buffer->Destroy();
      //
      pmrAllocatedBuffers.erase(it);          
    }    
  }

  void flush_pmr() {
  
    if (pmr_memory_pool) {
    
      pmrAllocatedBuffers.clear();
      pmrCachedBuffers.clear();      
      //
      for(auto &pmr_buffer : pmrBuffers) { 
        pmr_buffer->Destroy();  
        //        
        pmr_buffer.reset(nullptr);
      }
      //
      pmrBuffers.resize(0);
      
      pmr_memory_pool = false;
    }
  }
} // namespace pool


