#pragma once
#include <dslash.h>

template<typename Kernel, typename KernelArgs>
class Mat{
  private:
    std::unique_ptr<Kernel> dslash_kernel_ptr;

  public:

    Mat(const KernelArgs &args) : dslash_kernel_ptr(new Kernel(args)) {}

    void operator()(auto &out, auto &in){
      assert(in.GetFieldOrder() == FieldOrder::LexFieldOrder);
      // Take into account only internal points:
      const auto [Nx, Ny] = in.GetCBDims(); //Get CB dimensions

      auto X = std::views::iota(0, Nx-1);
      auto Y = std::views::iota(0, Ny-1);

      auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest
#if 0      
      //auto transformer = [](const auto &x, const auto a, const auto &y) {return (x-a*y);};
      auto DslashKernel = [&dslash_kernel = *dslash_kernel_ptr, transformer_ = transformer, out_ = out.Get(), in_ = in.Get()] (const auto i) { 
        //
        dslash_kernel.apply(transformer_, out_, in_, i); 
      };
#endif      

      auto DslashKernel = [&dslash_kernel = *dslash_kernel_ptr, out_ = out.Get(), in_ = in.Get()] (const auto i) { 
      
        const auto kappa = dslash_kernel.args.param.kappa;      
      
        auto transformer = [=](const auto &x, const auto &y) {return (x-kappa*y);};      
        //
        dslash_kernel.apply(transformer, out_, in_, i); 
      };
      //
      std::for_each(std::execution::par_unseq,
                    idx.begin(),
                    idx.end(),
                    DslashKernel);
    }
};


