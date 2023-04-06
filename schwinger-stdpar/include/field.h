#pragma once

#include <field_concepts.h>
#include <enums.h>
#include <assert.h>
#include <memory>

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
    //const std::array<int, type == FieldType::VectorFieldType ? nFace*ncolor*ncolor*ndim: nFace*nspin*ndim> comm_dirs;
    const std::array<int, nFace*ndim> comm_dir;

    const FieldOrder         order;        		
    const FieldSiteSubset    subset;
    const FieldParity        parity;

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
	            parity(parity){} 

    FieldDescriptor(const FieldDescriptor &args, const FieldSiteSubset subset,  const FieldParity parity) : 
	    dir{subset == FieldSiteSubset::ParitySiteSubset && args.subset == FieldSiteSubset::FullSiteSubset ? args.dir[0] / 2 : args.dir[0], args.dir[1]},
	    comm_dir{args.dir[1], subset == FieldSiteSubset::ParitySiteSubset && args.subset == FieldSiteSubset::FullSiteSubset ? args.dir[0] / 2 : args.dir[0]},
	    order(args.order),
	    subset(subset),
	    parity(parity){} 


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
    //template <GenericContainerTp T = container_tp, typename std::enable_if_t<!is_allocated_type_v<T>>* = nullptr>
    template <ContainerViewTp T>    
    Field(const T &src, const Arg &arg) : v(src), arg(arg) {}
    
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
    }
    //
    void destroy(){
      static_assert(is_allocated_type_v<container_tp>, "Cannot resize a non-owner field!");

      v.resize(0ul);
      ghost.resize(0ul);     
    }

    //Return a reference to the data container (adapter)
    auto& Data( ) { return v; }

    //Return a reference to the object (data access via container adapter )
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

