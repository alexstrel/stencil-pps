#pragma once

#include <assert.h>
#include <algorithm>
#include <ranges>
#include <execution>
#include <numeric>
#include <memory>
#include <memory_resource>

#include <fields/field_concepts.h>
#include <core/enums.h>
#include <core/memory.h>


template<std::size_t nD, std::size_t nS, std::size_t nC>
consteval FieldType get_field_type() {

  if constexpr (nD != invalid_dir and nS == invalid_spin  and nC != invalid_color){
    return FieldType::VectorFieldType;	  
  } else if constexpr (nD == invalid_dir and nS != invalid_spin  and nC != invalid_color) {
    return FieldType::SpinorFieldType;	  
  }

  return FieldType::InvalidFieldType;
}

template<std::size_t nDir = invalid_dir, std::size_t nSpin = invalid_spin, std::size_t nColor = invalid_color>
class FieldDescriptor {
  public: 
    static constexpr std::size_t ndim   = 2;
    static constexpr std::size_t ndir   = nDir;                    //vector field dim   (2 for U1 gauge)	  
    static constexpr std::size_t nspin  = nSpin;                   //number of spin dof (2 for spinor)
    static constexpr std::size_t ncolor = nColor;                  //for all fields

    static constexpr FieldType  type = get_field_type<ndir, nspin, ncolor>();

    static constexpr int nFace       = type == FieldType::VectorFieldType ? 1 : 2; 

    const std::array<int, ndim> dir;
    const std::array<int, nFace*ndim> comm_dir;

    const FieldOrder         order  = FieldOrder::InvalidFieldOrder;        		
    const FieldSiteSubset    subset = FieldSiteSubset::InvalidSiteSubset;
    const FieldParity        parity = FieldParity::InvalidFieldParity;

    std::shared_ptr<PMRBuffer> pmr_buffer;
    //
    bool is_exclusive;    
    //
    bool is_reserved;

    FieldDescriptor()                        = default;
    FieldDescriptor(const FieldDescriptor& ) = default;
    FieldDescriptor(FieldDescriptor&& )      = default;

    FieldDescriptor(const std::array<int, ndim> dir, 
                    const std::array<int, ndim*nFace> comm_dir,
	            const FieldOrder order         = FieldOrder::LexFieldOrder,
	            const FieldSiteSubset subset   = FieldSiteSubset::FullSiteSubset,  
	            const FieldParity parity       = FieldParity::InvalidFieldParity,
	            const bool is_exclusive        = true) : 
	            dir{dir},
                    comm_dir{comm_dir},
	            order(order),
	            subset(subset),
	            parity(parity), 
                    pmr_buffer(nullptr),
                    is_exclusive(is_exclusive),
                    is_reserved(false)  {} 

    FieldDescriptor(const FieldDescriptor &args, const FieldSiteSubset subset,  const FieldParity parity) : 
	            dir{subset == FieldSiteSubset::ParitySiteSubset && args.subset == FieldSiteSubset::FullSiteSubset ? args.dir[0] / 2 : args.dir[0], args.dir[1]},
	            comm_dir{args.dir[1], subset == FieldSiteSubset::ParitySiteSubset && args.subset == FieldSiteSubset::FullSiteSubset ? args.dir[0] / 2 : args.dir[0]},
	            order(args.order),
	            subset(subset),
	            parity(parity),
                    pmr_buffer(nullptr),
                    is_exclusive(args.IsExclusive()),
                    is_reserved(args.IsReserved())  {} 

    //Use it for block fields only:
    FieldDescriptor(const FieldDescriptor &args, 
                    const std::shared_ptr<PMRBuffer> extern_pmr_buffer) :
                    dir(args.dir),
                    comm_dir(args.comm_dir),
                    order(args.order),
                    subset(args.subset),
                    parity(args.parity),
                    pmr_buffer(extern_pmr_buffer),
                    is_exclusive(extern_pmr_buffer->IsExclusive()),
                    is_reserved(extern_pmr_buffer->IsReserved())  {}
       

    decltype(auto) GetFieldSize() const {
      if  constexpr (type == FieldType::ScalarFieldType) {
        return dir[0]*dir[1];
      } else if constexpr (type == FieldType::VectorFieldType) {
	return dir[0]*dir[1]*nDir*nColor*nColor;
      } else if constexpr (type == FieldType::SpinorFieldType) {
	return dir[0]*dir[1]*nSpin*nColor;
      }
      //
      return static_cast<std::size_t>(0);
    } 

    decltype(auto) GetGhostSize(int i, int face_idx = 0) const {
      if  constexpr (type == FieldType::ScalarFieldType) {
        return comm_dir[i*nFace+face_idx];
      } else if constexpr (type == FieldType::VectorFieldType) {
	return comm_dir[i*nFace+face_idx]*nColor*nColor;
      } else if constexpr (type == FieldType::SpinorFieldType) {
	return comm_dir[i*nFace+face_idx]*nSpin*nColor;
      }
      //
      return static_cast<std::size_t>(0);
    }

    decltype(auto) GetLatticeDims() const {
      return std::tie(dir[0], dir[1]);	    
    }

    decltype(auto) GetParityLatticeDims() const {
      const int xh = subset == FieldSiteSubset::FullSiteSubset ? dir[0] / 2 : dir[0];	    
      return std::make_tuple(xh, dir[1]);
    }
    
    template<ArithmeticTp T>
    void RegisterPMRBuffer(const bool is_reserved = false) {  
      // 
      const std::size_t nbytes = GetFieldSize()*sizeof(T);
      //
      if (pmr_buffer != nullptr) pmr_buffer.reset(); 
      //
      if ( is_exclusive ) {
        pmr_buffer = pmr_pool::pmr_malloc<true>(nbytes, is_reserved);
      } else {
        pmr_buffer = pmr_pool::pmr_malloc<false>(nbytes, is_reserved);      
      }

      this->is_reserved = is_reserved;
    }    

    void UnregisterPMRBuffer() {
      //
      if(pmr_buffer != nullptr) { 
        pmr_buffer->Release();
        //
        pmr_buffer.reset(); 
      }
    }
    
    void ResetPMRBuffer() {
      //
      if(pmr_buffer != nullptr) { 
        pmr_buffer.reset(); 
      }
    }

    void ReleasePMRBuffer() const { if(pmr_buffer != nullptr and not is_reserved) pmr_buffer->Release(); }     
    
    template<ArithmeticTp T>    
    bool IsReservedPMR() const {
      //
      const std::size_t nbytes = GetFieldSize()*sizeof(T);    
      //
      return pmr_buffer->IsReserved(nbytes);
    } 
    
    void SetExclusive()   { is_exclusive = true; }    
    void SetShared()      { is_exclusive = false; }        

    void SetReserved()    { is_reserved = true; }

    bool IsExclusive() const { return is_exclusive; }

    bool IsReserved()  const { return is_reserved; }

    auto GetState() const { 
      if (pmr_buffer != nullptr) {
        return pmr_buffer->State(); 
      }

      return PMRState::InvalidState;
    } 

    auto operator=(const FieldDescriptor&) -> FieldDescriptor& = default;
    auto operator=(FieldDescriptor&&     ) -> FieldDescriptor& = default;
};

template<int nDir,  int nColor = 1> using GaugeFieldArgs  = FieldDescriptor<nDir, invalid_spin, nColor>;

template<int nSpin, int nColor = 1> using SpinorFieldArgs = FieldDescriptor<invalid_dir,  nSpin, nColor>;

