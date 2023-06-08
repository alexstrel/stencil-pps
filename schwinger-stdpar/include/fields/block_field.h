#pragma once

#include <memory_resource>

template<GenericSpinorFieldTp spinor_tp, typename Arg, bool is_exclusive>
class BlockSpinor; // forward declare to make function definition possible

template <GenericSpinorFieldTp spinor_tp, typename Arg, bool use_pmr_buffer = false, bool is_exclusive = true>
decltype(auto) create_block_spinor(const Arg &arg_, const std::size_t n) {//offset for block spinor

  using data_tp = spinor_tp::container_tp::value_type;

  if constexpr ( use_pmr_buffer ) {
    const std::size_t pmr_bytes = (arg_.GetFieldSize()+arg_.GetGhostZoneSize())*sizeof(data_tp)*n;

    const bool reserved = true;
    
    auto pmr_buffer = pmr_pool::pmr_malloc<is_exclusive>(pmr_bytes, reserved);
    //
    auto pmr_arg = Arg{arg_, pmr_buffer};
    //
    return BlockSpinor<spinor_tp, Arg, is_exclusive>(pmr_arg, n, reserved);
  } else {
    auto arg = Arg{arg_};
  
    return BlockSpinor<spinor_tp, Arg>(arg, n);
  }
}

template<GenericSpinorFieldTp spinor_t, typename SpinorArg, bool is_exclusive = true>
class BlockSpinor{
  public:	
    using block_container_tp = std::vector<spinor_t> ;

    using container_tp  = typename spinor_t::container_tp;
    using arg_tp        = SpinorArg;
    using spinor_tp     = spinor_t;

    block_container_tp v;

    SpinorArg args;

    template<SpinorFieldTp T = spinor_t>
    BlockSpinor(const SpinorArg &args, const std::size_t n) : args(args) {
      v.reserve(n);
     
      for(int i = 0; i < n; i++) {
	v.push_back(create_field<container_tp, SpinorArg>(args));      
      }
    }
    
    template<PMRSpinorFieldTp T = spinor_t>
    BlockSpinor(const SpinorArg &args_, const std::size_t n, const bool is_reserved) : args(args_) {
      using data_tp = container_tp::value_type;;

      v.reserve(n);

      for(int i = 0; i < n; i++) {
        v.push_back(create_field_with_buffer<container_tp, SpinorArg, is_exclusive>(args, is_reserved));
      }
      //
      args.UpdatedReservedPMR();//now locked
    }

    template <SpinorFieldViewTp T = spinor_t>    
    BlockSpinor(const SpinorArg &args_, const std::size_t n) : args(args_) { v.reserve(n); }    

    decltype(auto) ConvertToView() {
      static_assert(is_allocator_aware_type<container_tp>, "Cannot reference a non-owner field!");
      
      using spinor_view_t = decltype(std::declval<spinor_t>().View());      

      const std::size_t n = nComponents();

      auto block_spinor_view = BlockSpinor<spinor_view_t, decltype(args), is_exclusive>{args, n};

      auto&& src_v = block_spinor_view.Get();

      for(auto &v_el : v) { src_v.push_back(v_el.View()); }

      return block_spinor_view;     
    }

    decltype(auto) ConvertToParityView(const FieldParity parity ) {
      static_assert(is_allocator_aware_type<container_tp>, "Cannot reference a non-owner field!");
      
      if (spinor_tp::nParity != 2) {
        std::cerr << "Cannot get a parity component from a non-full field, exiting...\n" << std::endl;
        std::quick_exit( EXIT_FAILURE );
      }
      
      using spinor_parity_view_t = decltype(std::declval<spinor_t>().ParityView(parity));

      const std::size_t n = nComponents();

      auto block_spinor_parity_view = BlockSpinor<spinor_parity_view_t, decltype(args), is_exclusive>{args, n};

      auto&& src_v = block_spinor_parity_view.Get();

      for(auto &v_el : v) { src_v.push_back(v_el.ParityView(parity)); }

      return block_spinor_parity_view;
    }

    auto ConvertToEvenView() { return ConvertToParityView(FieldParity::EvenFieldParity );}
    auto ConvertToOddView()  { return ConvertToParityView(FieldParity::OddFieldParity  );}
    
    auto EODecompose() {
      assert(spinor_tp::nParity == 2);

      return std::make_tuple(this->ConvertToEvenView(), this->ConvertToOddView());
    }    

    auto& Get()  { return v; }

    decltype(auto) BlockView() { return std::span{v}; }

    auto GetDims()       const { return args.GetLatticeDims(); }
    auto GetCBDims()     const { return args.GetParityLatticeDims(); }    

    auto GetFieldOrder()  const { return args.order; }    

    auto GetFieldSubset() const { return (spinor_tp::nParity == 2 ? FieldSiteSubset::FullSiteSubset : (spinor_tp::nParity == 1 ? FieldSiteSubset::ParitySiteSubset : FieldSiteSubset::InvalidSiteSubset)); }

    auto nComponents()   const { return v.size(); } 

    void destroy() {
      
      for(auto &spinor : v) spinor.destroy();
      
      args.ReleasePMRBuffer();
    }

    decltype(auto) ExportArg() { return args; }      

    spinor_t& operator[](const std::size_t i) { return v[i]; }
};



