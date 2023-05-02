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
#include <memory.h>


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

    FieldDescriptor()                        = default;
    FieldDescriptor(const FieldDescriptor& ) = default;
    FieldDescriptor(FieldDescriptor&& )      = default;

    FieldDescriptor(const std::array<int, ndim> dir, 
                    const std::array<int, ndim*nFace> comm_dir,
	            const FieldOrder order         = FieldOrder::LexFieldOrder,
	            const FieldSiteSubset subset   = FieldSiteSubset::FullSiteSubset,  
	            const FieldParity parity       = FieldParity::InvalidFieldParity,
	            const bool is_exclusive        = false) : 
	            dir{dir},
                    comm_dir{comm_dir},
	            order(order),
	            subset(subset),
	            parity(parity), 
                    pmr_buffer(nullptr),
                    is_exclusive(is_exclusive)  {} 

    FieldDescriptor(const FieldDescriptor &args, const FieldSiteSubset subset,  const FieldParity parity, const bool is_exclusive = false) : 
	            dir{subset == FieldSiteSubset::ParitySiteSubset && args.subset == FieldSiteSubset::FullSiteSubset ? args.dir[0] / 2 : args.dir[0], args.dir[1]},
	            comm_dir{args.dir[1], subset == FieldSiteSubset::ParitySiteSubset && args.subset == FieldSiteSubset::FullSiteSubset ? args.dir[0] / 2 : args.dir[0]},
	            order(args.order),
	            subset(subset),
	            parity(parity),
                    pmr_buffer(nullptr),
                    is_exclusive(is_exclusive)  {} 

    //Use it for block fields only:
    FieldDescriptor(const FieldDescriptor &args, 
                    const std::shared_ptr<PMRBuffer> extern_pmr_buffer) :
                    dir(args.dir),
                    comm_dir(args.comm_dir),
                    order(args.order),
                    subset(args.subset),
                    parity(args.parity),
                    pmr_buffer(extern_pmr_buffer),
                    is_exclusive(extern_pmr_buffer->IsExclusive())  {}
       

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
    void RegisterPMRBuffer(const std::size_t n = 1) {  
      // 
      const std::size_t bytes = GetFieldSize()*sizeof(T)*n;
      //
      bool is_reserved = (n > 1);  
      //
      auto new_pmr_buffer = pool::pmr_malloc<is_exclusive>(nbytes, is_reserved);
      
      if (pmr_buffer->State() == PMRState::Locked) {
      }    
      //
      pmr_buffer.reset(new_pmr_buffer);
    }    

    void UnregisterPMRBuffer() {
      //
      if(pmr_buffer != nullptr) { 
        pmr_buffer->Release();
        //
        pmr_buffer.reset(nullptr); 
      }
    }
    
    void ResetPMRBuffer() {
      //
      if(pmr_buffer != nullptr) { 
        pmr_buffer.reset(nullptr); 
      }
    }     
    
    bool IsReservedPMR(const std::size_t n = 1) const {
      //
      const std::size_t nbytes = GetFieldSize()*sizeof(T)*n;    
      //
      return pmr_buffer->IsReserved(nbytes);
    } 
    
    void SetExclusive() { is_exclusive = true; }    

    auto operator=(const FieldDescriptor&) -> FieldDescriptor& = default;
    auto operator=(FieldDescriptor&&     ) -> FieldDescriptor& = default;
};

template<int nDir,  int nColor = 1> using GaugeFieldArgs  = FieldDescriptor<nDir, invalid_spin, nColor>;

template<int nSpin, int nColor = 1> using SpinorFieldArgs = FieldDescriptor<invalid_dir,  nSpin, nColor>;

