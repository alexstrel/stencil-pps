#pragma once

#include <memory_resource>

template<SpinorFieldTp spinor_tp, typename Arg, bool is_exclusive>
class BlockSpinor; // forward declare to make function definition possible

template <SpinorFieldTp spinor_tp, typename Arg, bool use_pmr_buffer = false, bool is_exclusive = true>
decltype(auto) create_block_spinor(const Arg &arg_, const std::size_t n) {//offset for block spinor

  using data_tp = spinor_tp::container_tp::value_type;

  if constexpr ( use_pmr_buffer ) {
    const std::size_t pmr_bytes = arg_.GetFieldSize()*sizeof(data_tp)*n;

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


template<SpinorFieldTp spinor_t, typename SpinorArg, bool is_exclusive = true>
class BlockSpinor{
  public:
    using spinor_view_t = decltype(std::declval<spinor_t>().View());	 
    using container_tp  = typename spinor_t::container_tp;
    using arg_tp        = SpinorArg;

    std::vector<spinor_t> v;
    std::vector<spinor_view_t> w;

    SpinorArg args;

    template<SpinorFieldTp T = spinor_t>
    BlockSpinor(const SpinorArg &args, const std::size_t n) : args(args) {
      v.reserve(n);
      w.reserve(n);
     
      for(int i = 0; i < n; i++) {
	v.push_back(create_field<container_tp, SpinorArg>(args));      
	w.push_back(v[i].View());
      }
    }
    
    template<PMRSpinorFieldTp T = spinor_t>
    BlockSpinor(const SpinorArg &args_, const std::size_t n, const bool is_reserved) : args(args_) {
      using data_tp = container_tp::value_type;
      //constexpr bool is_reserved = true;

      v.reserve(n);
      w.reserve(n);

      for(int i = 0; i < n; i++) {
        v.push_back(create_field_with_buffer<container_tp, SpinorArg, is_exclusive>(args, is_reserved));
        //
        w.push_back(v[i].View());
      }
      //
      args.UpdatedReservedPMR();//now locked
    }

    auto View() { return std::span{w}; }

    auto GetDims() const { return args.GetLatticeDims(); }

    auto GetFieldOrder() const { return args.order; }    

    auto Size() const { return v.size(); } 

    void destroy() {
      
      for(auto &spinor : v) spinor.destroy();
      
      args.ReleasePMRBuffer();
    }

    decltype(auto) ExportArg() { return args; }      

    spinor_t& operator[](const std::size_t i) { return v[i]; }
};



