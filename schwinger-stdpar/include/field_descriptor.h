#pragma once

#include <field_concepts.h>
#include <enums.h>
#include <assert.h>
#include <memory>
#include <memory_resource>

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

    const FieldOrder         order;        		
    const FieldSiteSubset    subset;
    const FieldParity        parity;

    std::shared_ptr<std::byte[]>                          pmr_buffer;//needed only for pmr-used spinor
    std::size_t                                           pmr_bytes;
    std::shared_ptr<std::pmr::monotonic_buffer_resource>  pmr_pool; 
    //
    bool is_imported;    

    FieldDescriptor() = default;
    FieldDescriptor(const FieldDescriptor& ) = default;
    FieldDescriptor(FieldDescriptor&& )      = default;

    FieldDescriptor(const std::array<int, ndim> dir, 
                    const std::array<int, ndim*nFace> comm_dir,
	            const FieldOrder order         = FieldOrder::LexFieldOrder,
	            const FieldSiteSubset subset   = FieldSiteSubset::FullSiteSubset,  
	            const FieldParity parity       = FieldParity::InvalidFieldParity) : 
	            dir{dir},
                    comm_dir{comm_dir},
	            order(order),
	            subset(subset),
	            parity(parity), 
                    pmr_buffer(nullptr),
                    pmr_bytes(0ul), 
                    pmr_pool(nullptr), 
                    is_imported(false)  {} 

    FieldDescriptor(const FieldDescriptor &args, const FieldSiteSubset subset,  const FieldParity parity) : 
	            dir{subset == FieldSiteSubset::ParitySiteSubset && args.subset == FieldSiteSubset::FullSiteSubset ? args.dir[0] / 2 : args.dir[0], args.dir[1]},
	            comm_dir{args.dir[1], subset == FieldSiteSubset::ParitySiteSubset && args.subset == FieldSiteSubset::FullSiteSubset ? args.dir[0] / 2 : args.dir[0]},
	            order(args.order),
	            subset(subset),
	            parity(parity),
                    pmr_buffer(nullptr),
                    pmr_bytes(0ul), 
                    pmr_pool(nullptr),
                    is_imported(false)  {} 

    FieldDescriptor(const FieldDescriptor &args, 
                    const std::tuple<std::shared_ptr<std::byte[]>, std::size_t> extern_pmr) :
                    dir(args.dir),
                    comm_dir(args.comm_dir),
                    order(args.order),
                    subset(args.subset),
                    parity(args.parity),
                    pmr_buffer(std::get<0>(extern_pmr)),
                    pmr_bytes(std::get<1>(extern_pmr)),
                    pmr_pool(nullptr),
                    is_imported(true)  {}
       

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
    void AllocatePMRBuffer(const std::size_t n = 1) {
      if( is_imported ) return;

      const std::size_t bytes = GetFieldSize()*sizeof(T)*n;      

      if (pmr_buffer != nullptr && bytes <= pmr_bytes) return;
      
      if ( pmr_buffer != nullptr ) pmr_buffer.reset(); 
      
      pmr_bytes  = bytes;
      pmr_buffer = std::make_shared<std::byte[]>(pmr_bytes); 
    }

    template<ArithmeticTp T>
    void RegisterPMRPool(const std::size_t offset = 0ul) {
      if( pmr_buffer == nullptr) return;

      const std::size_t bytes = GetFieldSize()*sizeof(T);

      if(offset+bytes > pmr_bytes) return; 

      pmr_pool = std::make_shared< std::pmr::monotonic_buffer_resource >( pmr_buffer.get() + offset, bytes );
    }

    void ReleasePMR() {
      //
      if(pmr_pool   != nullptr) pmr_pool.reset();
      //
      if(pmr_buffer != nullptr and not is_imported) pmr_buffer.reset(); 

      pmr_bytes  = 0ul;
    }
 
    decltype(auto) ExportPMR() const { return std::tie(pmr_buffer, pmr_bytes); }

    void ImportPMR (const std::tuple<std::shared_ptr<std::byte[]>, std::size_t > &extern_pmr) {
      //
      if(pmr_pool != nullptr) {
        pmr_pool.reset();
      }

      if(pmr_buffer != nullptr and !is_imported) {
        pmr_buffer.reset();
        pmr_bytes  = 0ul;
      }
      
      pmr_buffer = std::get<std::shared_ptr<std::byte[]>>( extern_pmr );
      pmr_bytes  = std::get<std::size_t>( extern_pmr );
      //
      is_imported = true;
    }

    bool CheckPMRAllocation(const std::size_t bytes) const {
      return (pmr_buffer != nullptr and pmr_bytes >= bytes);
    }    

    auto operator=(const FieldDescriptor&) -> FieldDescriptor& = default;
    auto operator=(FieldDescriptor&&     ) -> FieldDescriptor& = default;
};

template<int nDir,  int nColor = 1> using GaugeFieldArgs  = FieldDescriptor<nDir, invalid_spin, nColor>;

template<int nSpin, int nColor = 1> using SpinorFieldArgs = FieldDescriptor<invalid_dir,  nSpin, nColor>;

