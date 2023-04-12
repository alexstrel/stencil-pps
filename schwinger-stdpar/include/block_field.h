#pragma once

template<SpinorFieldTp spinor_t, typename SpinorArg>
class BlockSpinor{
  public:
    using spinor_view_t = decltype(std::declval<spinor_t>().View());	 
    using container_tp  = typename spinor_t::container_tp;

    std::vector<spinor_t> v;
    std::vector<spinor_view_t> w;

    const SpinorArg args;

    BlockSpinor(const SpinorArg &args, const std::size_t n) : args(args) {
      v.reserve(n);
      w.reserve(n);
     
      for(int i = 0; i < n; i++) {
	v.push_back(create_field<container_tp, SpinorArg>(args));      
	w.push_back(v[i].View());
      }
    }
    
    template<PMRSpinorFieldTp T = spinor_t>
    BlockSpinor(const SpinorArg &args_, const std::size_t n) : args(args_) {
      using data_tp = container_tp::value_type;

      if( !std::is_same< typename T::container_tp::allocator_type , std::pmr::polymorphic_allocator<data_tp> >::value ) exit(-1);

      v.reserve(n);
      w.reserve(n);

      args.template AllocatePMRBuffer<data_tp>(n);       

      for(int i = 0; i < n; i++) {
        const std::size_t offset = i*args.GetFieldSize()*sizeof(data_tp);
        // 
        v.push_back(create_field_with_buffer<container_tp, SpinorArg>(args, offset));
        w.push_back(v[i].View());
      }
    }

    auto View() { return std::span{w}; }

    auto GetDims() const { return args.GetLatticeDims(); }

    auto GetFieldOrder() const { return args.order; }    

    auto Size() const { return v.size(); } 
};



