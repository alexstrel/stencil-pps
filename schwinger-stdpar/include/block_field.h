#pragma once

#include <memory_resource>

template<SpinorFieldTp spinor_tp, typename Arg>
class BlockSpinor; // forward declare to make function definition possible

template <PMRSpinorFieldTp pmr_spinor_tp, typename Arg>
decltype(auto) create_pmr_block_spinor(const Arg &arg_, const std::size_t n) {//offset for block spinor

  using data_tp = pmr_spinor_tp::container_tp::value_type;

  const std::size_t pmr_bytes = arg_.GetFieldSize()*sizeof(data_tp)*n;

  auto arg = Arg{arg_};

  if( arg.CheckPMRAllocation(pmr_bytes) == false ) {//just in case if the buffer is not allocated (or does not have an appropriate size)
    auto pmr_buffer = std::make_shared<std::byte[]>(pmr_bytes);
    arg.ImportPMR(std::tie(pmr_buffer, pmr_bytes));
  } 

  return BlockSpinor<pmr_spinor_tp, Arg>(arg, n);
}

template <PMRSpinorFieldTp pmr_spinor_t, typename Arg, PMRContainerTp pmr_container_tp = pmr_spinor_t::container_tp>
decltype(auto) export_pmr_block_spinor(BlockSpinor<pmr_spinor_t, Arg>& src_pmr_block_spinor, const std::size_t n, const bool reset_src = true) {//
  using data_tp = pmr_container_tp::value_type;
  //
  auto arg = Arg{src_pmr_block_spinor.ExportArg()};

  const std::size_t pmr_bytes = arg.GetFieldSize()*sizeof(data_tp)*n;
  
  if( arg.CheckPMRAllocation(pmr_bytes) == false ) {//just in case if the buffer is not allocated (or does not have an appropriate size)
    auto pmr_buffer = std::make_shared<std::byte[]>(pmr_bytes);
    arg.ImportPMR(std::tie(pmr_buffer, pmr_bytes));
  }  

  using new_pmr_spinor_tp = Field<pmr_container_tp, decltype(arg)>; 

  return BlockSpinor<new_pmr_spinor_tp, decltype(arg)>(arg, n); 
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

    void destroy() {
      for(auto &spinor : v) spinor.destroy();
    }

    decltype(auto) ExportArg() { return args; }      

    spinor_t& operator[](const std::size_t i) { return v[i]; }
};



