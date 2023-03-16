#pragma once

#include <common.h>
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
class FieldArgs {
  public: 
    static constexpr std::size_t ndir   = nDir;                    //vector field dim   (2 for U1 gauge)	  
    static constexpr std::size_t nspin  = nSpin;                   //number of spin dof (2 for spinor)
    static constexpr std::size_t ncolor = nColor;                  //for all fields

    static constexpr FieldType  type = get_field_type<ndir, nspin, ncolor>();

    const std::array<int, 2> dir;

    const FieldOrder         order;        		
    const FieldSiteSubset    subset;
    const FieldParity        parity;

    FieldArgs() = default;
    FieldArgs(const FieldArgs& ) = default;
    FieldArgs(FieldArgs&& )      = default;

    FieldArgs(const int L, 
	      const int T, 
	      const FieldOrder order         = FieldOrder::LexFieldOrder,
	      const FieldSiteSubset subset   = FieldSiteSubset::FullSiteSubset,  
	      const FieldParity parity       = FieldParity::InvalidFieldParity) : 
	      dir{L, T},
	      order(order),
	      subset(subset),
	      parity(parity){} 

    FieldArgs(const FieldArgs &args, const FieldSiteSubset subset,  const FieldParity parity) : 
	    dir{subset == FieldSiteSubset::ParitySiteSubset && args.subset == FieldSiteSubset::FullSiteSubset ? args.dir[0] / 2 : args.dir[0], args.dir[1]},
	    order(args.order),
	    subset(subset),
	    parity(parity){}  

    decltype(auto) GetFieldSize() const {
      if  constexpr (type == FieldType::ScalarFieldType) {
        return dir[0]*dir[1];
      } else if constexpr (type == FieldType::VectorFieldType) {
	return dir[0]*dir[1]*nDir*nColor;
      } else if constexpr (type == FieldType::SpinorFieldType) {
	return dir[0]*dir[1]*nSpin*nColor;
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

    auto operator=(const FieldArgs&) -> FieldArgs& = default;
    auto operator=(FieldArgs&&     ) -> FieldArgs& = default;
};

using GaugeFieldArgs  = FieldArgs<2, invalid_spin, 1>;
using SpinorFieldArgs = FieldArgs<invalid_dir,  2, 1>;

template<GenericContainerTp Ct, typename Arg>
class Field; // forward declare to make function definition possible
 
template <AllocatedContainerTp alloc_container_tp, typename Arg>
decltype(auto) create_field(const Arg &arg) {
  return Field<alloc_container_tp, Arg>(arg);
}

template <GenericContainerTp container_tp, typename Arg>
class Field{
  public:	
    //
    using data_tp = typename container_tp::value_type;    

    static constexpr std::size_t nDir    = Arg::ndir;
    static constexpr std::size_t nSpin   = Arg::nspin;                    
    static constexpr std::size_t nColor  = Arg::ncolor;                    

  private: 
    container_tp v;

    const Arg arg;//copy of the arguments

    Field(const Arg &arg) : v(arg.GetFieldSize()),
                            arg(arg){}

    template <AllocatedContainerTp alloc_container_tp, typename ArgTp>
    friend decltype(auto) create_field(const ArgTp &arg);

  public:

    Field()              = default;
    Field(const Field &) = default;
    Field(Field &&)      = default;    
    // 
    template <GenericContainerTp T = container_tp, typename std::enable_if_t<!is_allocated_type_v<T>>* = nullptr>
    Field(const T &src, const Arg &arg) : v(src), arg(arg) {}

    //Return a reference to the data container (adapter)
    auto& Data( ) { return v; }

    //Return a reference to the object (data access via container adapter )
    decltype(auto) Reference() {
      return Field<std::span<data_tp>, decltype(arg)>(std::span{v}, arg);	    
    }

    decltype(auto) ParityReference(const FieldParity parity ) {// return a reference to the parity component
      //
      if (arg.subset != FieldSiteSubset::FullSiteSubset) {
        std::cerr << "Cannot get a parity component from a non-full field, exiting...\n" << std::endl;
	exit(-1);
      }
      //
      auto parity_arg = FieldArgs(this->arg, FieldSiteSubset::ParitySiteSubset, parity);
      //
      const auto parity_length = GetParityLength();
      const auto parity_offset = parity == FieldParity::EvenFieldParity ? 0 : parity_length;

      return Field<std::span<data_tp>, decltype(parity_arg)>(std::span{v}.subspan(parity_offset, parity_length), parity_arg);
    }

    auto Even() { return ParityReference(FieldParity::EvenFieldParity );}
    auto Odd()  { return ParityReference(FieldParity::OddFieldParity  );}

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
    auto Accessor() const {
       //
       static_assert(nColor == 1, "Currently only O(1) model is supported.");

       using dyn_indx_type     = std::size_t;

       constexpr int nDoF = Arg::type == FieldType::VectorFieldType ? nDir*nColor*nColor : nSpin*nColor;

       using Dyn3DMap          = stdex::layout_stride::mapping<stdex::extents<dyn_indx_type, stdex::dynamic_extent, std::dynamic_extent, nDoF>>;
       using StridedDyn3DView  = stdex::mdspan<data_tp, stdex::extents<dyn_indx_type, stdex::dynamic_extent, stdex::dynamic_extent, nDoF>, stdex::layout_stride>;

       return StridedDyn3DView(v.data(), Dyn3DMap{stdex::extents<dyn_indx_type, stdex::dynamic_extent, stdex::dynamic_extent, nDoF>{arg.dir[0], arg.dir[1], nDoF}, std::array<dyn_indx_type, 3>{1, arg.dir[0], arg.dir[0]*arg.dir[1]}}) ;
    }
    
    //Direct field accessors (note that ncolor always 1, so no slicing for this dof):
    auto ExtAccessor() const {
       //
       static_assert(nColor == 1, "Currently only O(1) model is supported.");

       using dyn_indx_type     = std::size_t;

       constexpr int nDoF = Arg::type == FieldType::VectorFieldType ? nDir*nColor*nColor : nSpin*nColor;

       constexpr int nparity = 2;

       using Dyn3DMap          = stdex::layout_stride::mapping<stdex::extents<dyn_indx_type, stdex::dynamic_extent, std::dynamic_extent, nDoF, nparity>>;
       using StridedDyn3DView  = stdex::mdspan<data_tp, stdex::extents<dyn_indx_type, stdex::dynamic_extent, stdex::dynamic_extent, nDoF, nparity>, stdex::layout_stride>;

       return StridedDyn3DView(v.data(), Dyn3DMap{stdex::extents<dyn_indx_type, stdex::dynamic_extent, stdex::dynamic_extent, nDoF, nparity>{arg.dir[0], arg.dir[1], nDoF, npariry}, std::array<dyn_indx_type, 4>{1, arg.dir[0], arg.dir[0]*arg.dir[1], arg.dir[0]*arg.dir[1]*nDof}}) ;
    }    
    
};

/**
 *  Helper function that creates allocated containers
 */
#if 0
template <AllocatedContainerTp alloc_container_tp, typename Arg>
decltype(auto) make_field(const Arg &arg) {
  return Field<alloc_container_tp, Arg>(arg);
}
#endif
