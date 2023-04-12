#pragma once

#include <field_concepts.h>
#include <enums.h>
#include <assert.h>
#include <memory>
#include <memory_resource>
//
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

    const FieldOrder         order;        		
    const FieldSiteSubset    subset;
    const FieldParity        parity;

    std::shared_ptr<std::byte[]>                          pmr_buffer;//needed only for pmr-used spinor
    std::size_t                                           pmr_bytes;
    std::shared_ptr<std::pmr::monotonic_buffer_resource>  pmr_pool;     

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
                    pmr_pool(nullptr)  {} 

    FieldDescriptor(const FieldDescriptor &args, const FieldSiteSubset subset,  const FieldParity parity) : 
	            dir{subset == FieldSiteSubset::ParitySiteSubset && args.subset == FieldSiteSubset::FullSiteSubset ? args.dir[0] / 2 : args.dir[0], args.dir[1]},
	            comm_dir{args.dir[1], subset == FieldSiteSubset::ParitySiteSubset && args.subset == FieldSiteSubset::FullSiteSubset ? args.dir[0] / 2 : args.dir[0]},
	            order(args.order),
	            subset(subset),
	            parity(parity),
                    pmr_buffer(nullptr),
                    pmr_bytes(0ul), 
                    pmr_pool(nullptr)  {} 


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
      if( pmr_buffer != nullptr) return;//e.g., block spinors may need to call it multiple times

      pmr_bytes = GetFieldSize()*sizeof(T)*n;
      pmr_buffer = std::make_shared<std::byte[]>(pmr_bytes); 
    }

    template<ArithmeticTp T>
    void RegisterPMRBuffer(const std::size_t offset = 0ul) {
      if( pmr_buffer == nullptr) return;

      const std::size_t bytes = GetFieldSize()*sizeof(T);

      if(offset+bytes > pmr_bytes) return; 

      pmr_pool = std::make_shared< std::pmr::monotonic_buffer_resource >( pmr_buffer.get() + offset, bytes );
    }

    void ReleasePMRBuffer() const {
      //
      if(pmr_pool   != nullptr) pmr_pool.reset();
      pmr_pool = nullptr;
      //
      if(pmr_buffer != nullptr) pmr_buffer.reset();

      pmr_buffer = nullptr;
      pmr_bytes  = 0ul;
    }
    
    auto operator=(const FieldDescriptor&) -> FieldDescriptor& = default;
    auto operator=(FieldDescriptor&&     ) -> FieldDescriptor& = default;
};

template<int nDir,  int nColor = 1> using GaugeFieldArgs  = FieldDescriptor<nDir, invalid_spin, nColor>;
template<int nSpin, int nColor = 1> using SpinorFieldArgs = FieldDescriptor<invalid_dir,  nSpin, nColor>;

template<GenericContainerTp Ct, typename Arg>
class Field; // forward declare to make function definition possible
 
template <ContainerTp alloc_container_tp, typename Arg>
decltype(auto) create_field(const Arg &arg) {
  return Field<alloc_container_tp, Arg>(arg);
}

template <PMRContainerTp pmr_container_tp, typename Arg>
decltype(auto) create_field_with_buffer( const Arg &arg, const std::size_t offset = 0 ) {//offset for block spinors only

  using data_tp = pmr_container_tp::value_type;

  auto arg_ = Arg{arg};

  arg_.template AllocatePMRBuffer<data_tp>();
  arg_.template RegisterPMRBuffer<data_tp>(offset);

  auto& pmr_pool_ = *arg_.pmr_pool;

  return Field<pmr_container_tp, Arg>(pmr_pool_, arg_);
}

template <GenericContainerTp generic_container_tp, typename Arg>
class Field{
  public:	
    //
    using container_tp = generic_container_tp;        
    using data_tp      = typename container_tp::value_type;

    static constexpr std::size_t nDir    = Arg::ndir;
    static constexpr std::size_t nSpin   = Arg::nspin;                    
    static constexpr std::size_t nColor  = Arg::ncolor;                    

  private: 
    container_tp v;
    container_tp ghost;

    const Arg arg;//copy of the arguments

    Field(const Arg &arg) : v(arg.GetFieldSize()),
                            arg(arg){}

    template <ContainerTp alloc_container_tp, typename ArgTp>
    friend decltype(auto) create_field(const ArgTp &arg);

  public:

    Field()              = default;
    Field(const Field &) = default;
    Field(Field &&)      = default;    
    // 
    template <PMRContainerTp T = container_tp>
    Field(std::pmr::monotonic_buffer_resource &pmr_pool, const Arg &arg) : v(arg.GetFieldSize(), &pmr_pool), arg(arg) { }
    //
    template <ContainerViewTp T>    
    Field(const T &src, const Arg &arg) : v{src}, arg(arg) {}
    
    // Needed for block-la operations
    constexpr std::size_t size() const { return 1ul; } 
    
    // Needed for block-la operations
    decltype(auto) operator[](int i) const { 
      if(i != 0) exit(-1); 
      //
      return *this; 
    }        
    //
    void move() {
//printf("PTR: %p, %d\n", this->arg.pmr_buffer.get(), this->arg.pmr_buffer.use_count());
//arg.ReleasePMRBuffer();
//printf("PTR: %p, %d\n", this->arg.pmr_buffer.get(), this->arg.pmr_buffer.use_count());
    }
    //
    void destroy(){
      static_assert(is_allocated_type_v<container_tp>, "Cannot resize a non-owner field!");

      v.resize(0ul);
      ghost.resize(0ul);     
    }

    //Return a reference to the data container
    auto& Data( ) { return v; }

    //Return a reference to the object (data access via std::span )
    decltype(auto) View() {
      static_assert(is_allocated_type_v<container_tp>, "Cannot reference a non-owner field!");

      return Field<std::span<data_tp>, decltype(arg)>(std::span{v}, arg);	    
    }

    decltype(auto) ParityView(const FieldParity parity ) {// return a reference to the parity component
      static_assert(is_allocated_type_v<container_tp>, "Cannot reference a non-owner field!");
      //
      if (arg.subset != FieldSiteSubset::FullSiteSubset) {
        std::cerr << "Cannot get a parity component from a non-full field, exiting...\n" << std::endl;
	exit(-1);
      }
      //
      auto parity_arg = FieldDescriptor(this->arg, FieldSiteSubset::ParitySiteSubset, parity);
      //
      const auto parity_length = GetParityLength();
      const auto parity_offset = parity == FieldParity::EvenFieldParity ? 0 : parity_length;

      return Field<std::span<data_tp>, decltype(parity_arg)>(std::span{v}.subspan(parity_offset, parity_length), parity_arg);
    }

    auto Even() { return ParityView(FieldParity::EvenFieldParity );}
    auto Odd()  { return ParityView(FieldParity::OddFieldParity  );}

    auto EODecompose() {
      assert(arg.subset == FieldSiteSubset::FullSiteSubset);	    

      return std::make_tuple(this->Even(), this->Odd());
    }

    auto GetLength()       const { return v.size(); }
    auto GetBytes()        const { return v.size()*sizeof(data_tp); }
    auto GetParityLength() const { return v.size() / (arg.subset == FieldSiteSubset::FullSiteSubset ? 2 : 1); }
    auto GetDims()         const { return arg.GetLatticeDims(); }
    auto GetCBDims()       const { return arg.GetParityLatticeDims(); }
   
    auto GetFieldSubset()  const { return arg.subset; }
    auto GetFieldOrder()   const { return arg.order; }    

    //Direct field accessors (note that ncolor always 1, so no slicing for this dof):
    template<bool is_constant = false>
    auto Accessor() const {
       //
       static_assert(nColor == 1, "Currently only O(1) model is supported.");

       using dyn_indx_type     = std::size_t;

       constexpr dyn_indx_type nDoF = Arg::type == FieldType::VectorFieldType ? nDir*nColor*nColor : nSpin*nColor;

       using Dyn3DMap  = stdex::layout_stride::mapping<stdex::extents<dyn_indx_type, stdex::dynamic_extent, std::dynamic_extent, nDoF>>;
       using Extents3D = stdex::extents<dyn_indx_type, stdex::dynamic_extent, stdex::dynamic_extent, nDoF>;
       using Strides3D = std::array<dyn_indx_type, 3>;       
       
       if constexpr (is_constant){
         return stdex::mdspan<const data_tp, Extents3D, stdex::layout_stride>{
                    v.data(), Dyn3DMap{Extents3D{arg.dir[0], arg.dir[1]}, Strides3D{1, arg.dir[0], arg.dir[0]*arg.dir[1]}}} ;
       } else {
         return stdex::mdspan<data_tp, Extents3D, stdex::layout_stride>{
                   v.data(), Dyn3DMap{Extents3D{arg.dir[0], arg.dir[1]}, Strides3D{1, arg.dir[0], arg.dir[0]*arg.dir[1]}}};
       }
    }
    
    //Direct field accessors (note that ncolor always 1, so no slicing for this dof):
    template<bool is_constant = false>    
    auto ExtAccessor() const {
       //
       static_assert(nColor == 1, "Currently only O(1) model is supported.");

       using dyn_indx_type     = std::size_t;

       constexpr dyn_indx_type nDoF = Arg::type == FieldType::VectorFieldType ? nDir*nColor*nColor : nSpin*nColor;

       constexpr dyn_indx_type nparity = 2;

       using Dyn4DMap  = stdex::layout_stride::mapping<stdex::extents<dyn_indx_type, stdex::dynamic_extent, std::dynamic_extent, nDoF, nparity>>;
       using Extents4D = stdex::extents<dyn_indx_type, stdex::dynamic_extent, stdex::dynamic_extent, nDoF, nparity>;
       using Strides4D = std::array<dyn_indx_type, 4>;       
       
       if constexpr (is_constant){
         return stdex::mdspan<const data_tp, Extents4D, stdex::layout_stride>{
                   v.data(), Dyn4DMap{Extents4D{arg.dir[0], arg.dir[1]}, Strides4D{1, arg.dir[0], arg.dir[0]*arg.dir[1], arg.dir[0]*arg.dir[1]*nDoF}}} ;
       } else {
         return stdex::mdspan<data_tp, Extents4D, stdex::layout_stride>{
                   v.data(), Dyn4DMap{Extents4D{arg.dir[0], arg.dir[1]}, Strides4D{1, arg.dir[0], arg.dir[0]*arg.dir[1], arg.dir[0]*arg.dir[1]*nDoF}}} ;       
       }
    }    
    
};

