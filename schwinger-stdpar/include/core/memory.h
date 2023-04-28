#pragma once
#include <common.h>
#include <memory>
#include <memory_resource>

static std::pmr::memory_resource *upstream = std::pmr::get_default_resource(); 

class PMRBuffer{
  public:
    //
    std::shared_ptr<std::byte[]>    pmr_buffer;//needed only for pmr-used spinor
    std::size_t                     pmr_bytes;

    PMRState  pmr_state = PMRState::InvalidState; 

    PMRBuffer(const std::size_t nbytes, const enum state) pmr_bytes(nbytes), pmr_state(state) : {
      pmr_buffer = std::make_shared<std::byte[]>(pmr_bytes);
    } 

    PMRBuffer()                             = default;
    PMRBuffer(const PMRBuffer &)            = default;
    PMRBuffer(PMRBuffer &&)                 = default;

    virtual ~PMRBuffer()                    = default;

    auto Data()  const { return pmr_buffer;       }
    auto Get()   const { return pmr_buffer.get(); }
    auto Byte()  const { return pmr_bytes;        }    
    auto State() const { return pmr_state;        }        
    
    auto SetState(PMRState state) { rpmr_state = state; }            

    PMRBuffer &operator=(const PMRBuffer &) = default;
    PMRBuffer &operator=(PMRBuffer &&)      = default;
};



inline void allocate_pmr_pool(const std::size_t bytes) {
  //
  pool_bytes = bytes;
  //
  if(pool_ptr == nullptr) pool_ptr = std::make_shared<std::byte[]>(pool_bytes);  
} 

inline decltype(auto) allocate_extern_pmr_pool(const std::size_t bytes) {
  return std::make_shared<std::byte[]>(bytes);
}
        
inline void release_pmr_pool() {
  if(pool_ptr != nullptr) pool_ptr.reset();

  pool_ptr  = nullptr;
  pool_bytes= 0ul;   
}

inline void release_extern_pmr_pool(std::shared_ptr<std::byte[]> &ptr) {
  if(ptr != nullptr) ptr.reset();
  ptr  = nullptr;
}

decltype(auto) get_pmr_pool(const size_t bytes, const std::size_t offset = 0ul) {
  if(pool_ptr == nullptr) {
    std::cerr << "Cannot create a pmr resource, upstream pool does not exist." << std::endl;
    exit(-1);
  }

  if (bytes > pool_bytes) {
    std::cerr << "Cannot create a pmr resource, requested number of pool bytes to big." << std::endl;
    exit(-1); 
  }
 
  return std::pmr::monotonic_buffer_resource{pool_ptr.get() + offset, bytes};
}

template<PMRContainerTp T>
decltype(auto) get_pmr_container(const size_t n) {
  if(pool_ptr == nullptr) {
    std::cerr << "Cannot create a pmr resource, upstream pool does not exist." << std::endl;
    exit(-1);
  }

  const std::size_t bytes = n*sizeof(typename T::value_type);

  if (bytes > pool_bytes) {
    std::cerr << "Cannot create a pmr resource, requested number of pool bytes to big." << std::endl;
    exit(-1);
  }

  auto pmr_pool = std::pmr::monotonic_buffer_resource{pool_ptr.get(), bytes};

  return T(n, &pmr_pool);
}

/**
    Derived from the QUDA library, original developer : Kate Clark (NVIDIA)
*/



namespace pmr_pool
{
  /** Cache of pmr-memory allocations.  We cache pmr
      memory allocations so that fields can reuse them.*/
  static std::multimap<std::size_t, PMRBuffer> pmrBuffers;

  /** Sizes of active pmr-memory allocations.  For convenience,
      we keep track of the sizes of active allocations (i.e., those not
      in the cache). */
  static std::map<std::byte*, PMRState> pmrReservedBuffers;

  static bool pmr_pool_init   = false;

  /** whether to use a memory pool allocator for pinned memory */
  static bool pmr_memory_pool = false;

  //void init() { pmr_pool_init = true; }

  template<bool is_exclusive = true>
  decltype(auto) pmr_malloc_(const std::size_t nbytes, const bool is_reserved = false)
  {
    //
    constexpr PRMState buffer_state = is_exclusive ? PMRState::Locked : PMRState::NonVacant;
    
    if ( pmrBuffers.empty() ) {
      //
      auto pmr_buffer = PMRBuffer{nbytes, is_reserved ? PMRState::Reserved : buffer_state};
      //
      pmrBuffers.insert(std::make_pair(nbytes, pmr_buffer));	
      //
      if ( is_reserved ) pmrReservedBuffers.insert(std::make_pair(pmr_buffer.Get(), buffer_state));
        
      return std::pmr::monotonic_buffer_resource{pmr_buffer.Get(), nbytes};
        
    } else {

      auto range_begin = pmrBuffers.lower_bound(nbytes);
      auto range_end   = pmrBuffers.end();
          
      if ( range_begin != range_end ) {
        //
        auto subrange_view = std::range::subrange(range_begin, range_end);
	  
        for(auto& [bytes, pmr_buffer] : subrange_view) {
	  if (pmr_buffer.State() == PMRState::Vacant) {
	     
	    pmr_buffer.SetState( is_reserved ? PMRState::Reserved : buffer_state );
	    //
            if ( is_reserved ) pmrReservedBuffers.insert(std::make_pair(pmr_buffer.Get(), buffer_state));	      
	    //
	    return std::pmr::monotonic_buffer_resource{pmr_buffer.Get(), nbytes};
	  } 
        } 
      }
    }//
    
    return std::pmr::monotonic_buffer_resource{nbytes, upstream};//in all other cases?
  }    

    void pinned_free_(const char *func, const char *file, int line, void *ptr)
    {
      if (pinned_memory_pool) {
        if (!pinnedSize.count(ptr)) { errorQuda("Attempt to free invalid pointer"); }
        pinnedCache.insert(std::make_pair(pinnedSize[ptr], ptr));
        pinnedSize.erase(ptr);
      } else {
        quda::host_free_(func, file, line, ptr);
      }
    }

    void *device_malloc_(const char *func, const char *file, int line, size_t nbytes)
    {
      void *ptr = nullptr;
      if (device_memory_pool) {
        if (deviceCache.empty()) {
          ptr = quda::device_malloc_(func, file, line, nbytes);
        } else {
          auto it = deviceCache.lower_bound(nbytes);
          if (it != deviceCache.end()) { // sufficiently large allocation found
            nbytes = it->first;
            ptr = it->second;
            deviceCache.erase(it);
          } else { // sacrifice the smallest cached allocation
            it = deviceCache.begin();
            ptr = it->second;
            deviceCache.erase(it);
            quda::device_free_(func, file, line, ptr);
            ptr = quda::device_malloc_(func, file, line, nbytes);
          }
        }
        deviceSize[ptr] = nbytes;
      } else {
        ptr = quda::device_malloc_(func, file, line, nbytes);
      }
      return ptr;
    }

    void device_free_(const char *func, const char *file, int line, void *ptr)
    {
      if (device_memory_pool) {
        if (!deviceSize.count(ptr)) { errorQuda("Attempt to free invalid pointer"); }
        deviceCache.insert(std::make_pair(deviceSize[ptr], ptr));
        deviceSize.erase(ptr);
      } else {
        quda::device_free_(func, file, line, ptr);
      }
    }

#ifdef NVSHMEM_COMMS
    void *shmem_malloc_(const char *func, const char *file, int line, size_t nbytes)
    {
      return quda::shmem_malloc_(func, file, line, nbytes);
    }

    void shmem_free_(const char *func, const char *file, int line, void *ptr)
    {
      quda::shmem_free_(func, file, line, ptr);
    }
#endif

    void flush_pinned()
    {
      if (pinned_memory_pool) {
        for (auto it : pinnedCache) { host_free(it.second); }
        pinnedCache.clear();
      }
    }

    void flush_device()
    {
      if (device_memory_pool) {
        for (auto it : deviceCache) { device_free(it.second); }
        deviceCache.clear();
      }
    }

  } // namespace pool




