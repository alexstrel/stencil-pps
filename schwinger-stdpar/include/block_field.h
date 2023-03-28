#pragma once

template<SpinorFieldTp spinor_t, typename SpinorArg>
class BlockSpinor{
  public:
    using spinor_view_t = decltype(std::declval<spinor_t>().View());	 
    using container_tp  = typename spinor_t::container_tp;

    std::vector<spinor_t> v;
    std::vector<spinor_view_t> w;

    const SpinorArg args;

    //template<GenericContainerTp container_tp>
    BlockSpinor(const SpinorArg &args, const std::size_t n) : args(args) {
      v.reserve(n);
      w.reserve(n);
      //
      for(int i = 0; i < n; i++) {
	v.push_back(create_field<container_tp, SpinorArg>(args));      
	w.push_back(v[i].View());
      }
    }

    auto View() { return std::span{w}; }

    auto GetDims() const { return args.GetLatticeDims(); }

    auto GetFieldOrder() const { return args.order; }    
};



