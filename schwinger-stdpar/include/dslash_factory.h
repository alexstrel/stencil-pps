#pragma once
#include <dslash.h>
#include <optional>

template<typename KernelArgs, template <typename Args> class Kernel>
class Mat{
  private:
    std::unique_ptr<Kernel<KernelArgs>> dslash_kernel_ptr;
    
  public:

    Mat(const KernelArgs &args) : dslash_kernel_ptr(new Kernel<KernelArgs>(args)) {}
   
    void operator()(auto &out, auto &in, auto&& transformer = std::nullopt){
      assert(in.GetFieldOrder() == FieldOrder::LexFieldOrder);
      // Extract dims:
      const auto [Nx, Ny] = in.GetCBDims(); //Get CB dimensions

      auto X = std::views::iota(0, Nx-1);
      auto Y = std::views::iota(0, Ny-1);

      auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest

      auto DslashKernel = [&dslash_kernel   = *dslash_kernel_ptr, 
	                   transformer      = transformer, 
			   out_             = out.Get(), 
			   in_              = in.Get()           ] (const auto coords) { 
                             //
                             dslash_kernel.apply(transformer, out_, in_, coords); 
                           };
      //
      std::for_each(std::execution::par_unseq,
                    idx.begin(),
                    idx.end(),
                    DslashKernel);
    }
};


