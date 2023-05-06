#pragma once

#include <fields/field_descriptor.h>

template<GenericContainerTp Ct, typename Arg>
class Field; // forward declare to make function definition possible
 
template <ContainerTp alloc_container_tp, typename Arg>
decltype(auto) create_field(const Arg &arg) {
  return Field<alloc_container_tp, Arg>(arg);
}

template <FieldTp field_tp, ContainerTp dst_container_tp = field_tp::container_tp, bool do_copy = false>
decltype(auto) create_field(field_tp &src) {
  //
  using src_container_tp = field_tp::container_tp;
  
  using src_data_tp      = field_tp::data_tp;  
  using dst_data_tp      = dst_container_tp::value_type;    
  //
  auto arg = src.ExportArg();
  //
  auto dst = Field<dst_container_tp, decltype(arg)>(arg);
  
  if constexpr (do_copy){
    //
    auto &&dst_view = dst.View();
    auto &&src_view = src.View();
    //    
    if constexpr (std::is_same_v< src_container_tp, dst_container_tp >) {
      std::copy( std::execution::par_unseq, dst_view.Begin(), dst_view.End(), src_view.Begin());
    } else {
      std::transform(std::execution::par_unseq, src_view.Begin(), src_view.End(), dst_view.Begin(), [=](const auto &in) { return static_cast<dst_data_tp>(in); } );
    }
  }
  //
  return dst;
}

template <PMRContainerTp pmr_container_tp, typename Arg>
decltype(auto) create_field_with_buffer(const Arg &arg_, const bool use_reserved = false) {
  using data_tp = pmr_container_tp::value_type;

  auto arg = Arg{arg_};

  if ( not use_reserved) {// if not use a reserved buffer, register new one
    arg.template RegisterPMRBuffer<data_tp>();
  } else {
  // use reserved, e.g., by a block field constructor, 
  // arg must have a valid pmr_buffer with PMRStatus::Reserved
    if ( not arg.template IsReservedPMR<data_tp>() ) {
      std::cerr << "Incorrect PMR buffer state, check reservation." << std::endl;
      exit(-1);    
    }
    arg.SetReserved();
  }
  
  if ( not arg.IsExclusive() ) std::cout << "Warning: creating non-exclusive PMR field." << std::endl;

  auto& pmr_pool_handle = *arg.pmr_buffer->Pool();

  return Field<pmr_container_tp, Arg>(pmr_pool_handle, arg);
}

template <GenericContainerTp generic_container_tp, typename Arg>
class Field{
  public:	
    //
    using container_tp  = generic_container_tp;        
    using data_tp       = typename container_tp::value_type;

    static constexpr std::size_t nDir    = Arg::ndir;
    static constexpr std::size_t nSpin   = Arg::nspin;                    
    static constexpr std::size_t nColor  = Arg::ncolor;                    

  private: 
    container_tp v;
    container_tp ghost;

    const Arg arg;//copy of the arguments

    Field(const Arg &arg) : v(arg.GetFieldSize()),
                            arg(arg){ }
    // 
    template <PMRContainerTp T = container_tp>
    explicit Field(std::pmr::monotonic_buffer_resource &pmr_pool, const Arg &arg) : v(arg.GetFieldSize(), &pmr_pool), arg(arg) { }

    template <ContainerTp alloc_container_tp, typename ArgTp>
    friend decltype(auto) create_field(const ArgTp &arg);    
    //    
    template <FieldTp field_tp, ContainerTp container_tp, bool do_copy>
    friend decltype(auto) create_field(field_tp &src);

    template <PMRContainerTp pmr_container_tp, typename ArgTp>
    friend decltype(auto) create_field_with_buffer(const ArgTp &arg_, const bool use_reserved);

  public:

    Field()              = default;
    Field(const Field &) = default;
    Field(Field &&)      = default;    
    // 
    template <ContainerViewTp T>    
    explicit Field(const T &src, const Arg &arg) : v{src}, arg(arg) {}
    
    // Needed for block-la operations
    constexpr std::size_t size() const { return 1ul; } 
    
    // Needed for block-la operations
    decltype(auto) operator[](int i) const { 
      if(i != 0) exit(-1); 
      //
      return *this; 
    }        
    //
    void show() {
      printf("PMR buffer pointer: %p, PMR buffer use count %d\n", this->arg.pmr_buffer.get(), this->arg.pmr_buffer.use_count());
    }
    
    void destroy() {
      static_assert(is_allocated_type_v<container_tp>, "Cannot resize a non-owner field!");

      v.resize(0ul);
      ghost.resize(0ul); 
      //
      if constexpr ( std::is_same_v< typename container_tp::allocator_type,  std::pmr::polymorphic_allocator<data_tp> > ) {
        arg.ReleasePMRBuffer();
      }
    }

    decltype(auto) ExportArg() { return Arg{arg}; }

    auto Begin() { return v.begin(); }
    auto End()   { return v.end();   }

    //Return a reference to the data container
    auto& Data( ) { return v; }

    //Return the data raw ptr
    auto Get( ) const { return v.data(); }    

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
    
    
    template<bool is_constant = false>
    auto FlatAccessor(const std::size_t stride = 1) const {
       //
       static_assert(nColor == 1, "Currently only O(1) model is supported.");
    
       using dyn_indx_type = std::size_t;

       using Dyn1DMap  = stdex::layout_stride::mapping<stdex::extents<dyn_indx_type, stdex::dynamic_extent>>;
       using Extents1D = stdex::extents<dyn_indx_type, stdex::dynamic_extent>;
       using Strides1D = std::array<dyn_indx_type, 1>;       
       
       const dyn_indx_type len = GetLength() / stride;
       
       if constexpr (is_constant){
         return stdex::mdspan<const data_tp, stdex::extents<dyn_indx_type, stdex::dynamic_extent>, stdex::layout_stride>{
                    v.data(), Dyn1DMap{Extents1D{len}, Strides1D{stride}}} ;
       } else {
         return stdex::mdspan<data_tp, stdex::extents<dyn_indx_type, stdex::dynamic_extent>, stdex::layout_stride>{
                    v.data(), Dyn1DMap{Extents1D{len}, Strides1D{stride}}} ;       
       }    
    }    
    
};

