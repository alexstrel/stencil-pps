#pragma once

#include <common.h>
#include <enums.h>

template<int nDir = invalid_dir, int nSpin = invalid_spin, int nColor = invalid_color>
class FieldArgs {
  public: 
    static constexpr int ndir   = nDir;                    //vector field dim   (2 for U1 gauge)	  
    static constexpr int nspin  = nSpin;                   //number of spin dof (2 for spinor)
    static constexpr int ncolor = nColor;                  //for all fields

    static constexpr FieldType  type = nDir > 1 ? FieldType::VectorFieldType : FieldType::SpinorFieldType;

    const std::array<int, 2> dir;		
    const FieldSiteSubset    subset;
    //
    const FieldParity        parity;

    FieldArgs() = default;
    FieldArgs(const FieldArgs& ) = default;
    FieldArgs(FieldArgs&& )      = default;

    FieldArgs(const int L, 
	      const int T, 
	      const FieldSiteSubset subset   = FieldSiteSubset::FullSiteSubset,  
	      const FieldParity parity       = FieldParity::InvalidFieldParity) : 
	      dir{L, T},
	      subset(subset),
	      parity(parity){} 

    FieldArgs(const FieldArgs &args, const FieldSiteSubset subset,  const FieldParity parity) : 
	    dir{subset == FieldSiteSubset::ParitySiteSubset && args.subset == FieldSiteSubset::FullSiteSubset ? args.dir[0] / 2 : args.dir[0], args.dir[1]},
	    subset(subset),
	    parity(parity){}  

    auto GetFieldSize() const {
      if  constexpr (type == FieldType::ScalarFieldType) {
        return dir[0]*dir[1];
      } else if constexpr (type == FieldType::VectorFieldType) {
	return dir[0]*dir[1]*nDir*nColor;
      } else if constexpr (type == FieldType::SpinorFieldType) {
	return dir[0]*dir[1]*nSpin*nColor;
      }
      //
      return 0;
    }

    auto operator=(const FieldArgs&) -> FieldArgs& = default;
    auto operator=(FieldArgs&&     ) -> FieldArgs& = default;
};

using GaugeFieldArgs  = FieldArgs<2, invalid_spin, 1>;
using SpinorFieldArgs = FieldArgs<invalid_dir,  2, 1>;

template <GenericContainerTp container_tp, typename Arg>
class Field{
  public:	
    //
    static constexpr int nDir    = Arg::ndir;
    static constexpr int nSpin   = Arg::nspin;                    
    static constexpr int nColor  = Arg::ncolor;                    

  private: 
    container_tp v;

    const Arg arg;//copy of the arguments

  public:
    //
    Field(const Arg &arg) : v(arg.GetFieldSize()), 
	                    arg(arg) {}
    Field(const container_tp &src, const Arg &arg) : v(src),
                            arg(arg) {}

    auto& Get( ) { return v; }

    auto Even() {
      if (arg.subset != FieldSiteSubset::FullSiteSubset) {
        std::cerr << "Cannot get a parity component from a non-full field\n" << std::endl;
      }
      using data_tp = typename container_tp::value_type;
      //
      auto even_arg = FieldArgs(this->arg, FieldSiteSubset::ParitySiteSubset, FieldParity::EvenFieldParity);
      //
      return Field<std::span<data_tp>, decltype(even_arg)>(std::span{v}.subspan(0, GetParityLength()), even_arg);
    }

    auto Odd() {
      if (arg.subset != FieldSiteSubset::FullSiteSubset) {
        std::cerr << "Cannot get a parity component from a non-full field\n" << std::endl;
      }
      using data_tp = typename container_tp::value_type;
      //
      auto odd_arg = FieldArgs(this->arg, FieldSiteSubset::ParitySiteSubset, FieldParity::OddFieldParity);
      //
      return Field<std::span<data_tp>, decltype(odd_arg)>(std::span{v}.subspan(GetParityLength(), GetParityLength()), odd_arg);
    }

    auto GetLength()       const { return v.size(); }
    auto GetParityLength() const { return v.size() / (arg.subset == FieldSiteSubset::FullSiteSubset ? 2 : 1); }

    //Field accessors (note that ncolor always 1, so no slicing for this dof):
    auto Accessor() {
       //
       static_assert(nColor == 1, "Currently only O(1) model is supported\n");

       using data_tp           = typename container_tp::value_type;
       using dyn_indx_type     = int;

       constexpr int nDoF = Arg::type == FieldType::VectorFieldType ? nDir*nColor : nSpin*nColor;

       using Dyn3DMap          = stdex::layout_stride::mapping<stdex::extents<dyn_indx_type, stdex::dynamic_extent, std::dynamic_extent, nDoF>>;
       using StridedDyn3DView  = stdex::mdspan<data_tp, stdex::extents<dyn_indx_type, stdex::dynamic_extent, stdex::dynamic_extent, nDoF>, stdex::layout_stride>;

       return StridedDyn3DView(v.data(), Dyn3DMap{stdex::extents<dyn_indx_type, stdex::dynamic_extent, stdex::dynamic_extent, nDir>{arg.dir[0], arg.dir[1], nDoF}, std::array<dyn_indx_type, 3>{1, arg.dir[0], arg.dir[0]*arg.dir[1]}}) ;
    }
};

