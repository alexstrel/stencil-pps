#pragma once

#include <memory_resource>

template<SpinorFieldTp spinor_tp, typename Arg>
class BlockSpinor; // forward declare to make function definition possible

template <PMRSpinorFieldTp pmr_spinor_tp, typename Arg>
decltype(auto) create_block_spinor_with_buffer(const Arg &arg_, const std::size_t n) {//offset for block spinor

  using data_tp = pmr_spinor_tp::container_tp::value_type;

  const std::size_t new_pmr_bytes = arg_.GetFieldSize()*sizeof(data_tp)*n;

  auto arg = Arg{arg_};

  if( arg.CheckPMR(new_pmr_bytes) == false ) {
printf("Realloc new buffer %d\n", new_pmr_bytes);
    auto new_pmr_buffer = std::make_shared<std::byte[]>(new_pmr_bytes);
    arg.ImportPMR(std::tie(new_pmr_buffer, new_pmr_bytes));
  } 

  //auto arg = Arg{arg_, std::tie(pmr_ptr, pmr_bytes)}; 
  //auto arg = Arg{arg_};
  //arg.template AllocatePMRBuffer<data_tp>(n);

  return BlockSpinor<pmr_spinor_tp, Arg>(arg, n);
}


template<SpinorFieldTp spinor_t, typename SpinorArg>
class BlockSpinor{
  public:
    using spinor_view_t = decltype(std::declval<spinor_t>().View());	 
    using container_tp  = typename spinor_t::container_tp;

    std::vector<spinor_t> v;
    std::vector<spinor_view_t> w;

    const SpinorArg args;

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
    BlockSpinor(const SpinorArg &args_, const std::size_t n) : args(args_) {
      using data_tp = container_tp::value_type;

      v.reserve(n);
      w.reserve(n);

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

    

    spinor_t& operator[](const std::size_t i) { return v[i]; }
};



