#pragma once
#include <dslash.h>
#include <optional>

template<typename KernelArgs, template <typename Args> class Kernel, typename TransformParams>
class Mat{
  private:
    std::unique_ptr<Kernel<KernelArgs>> dslash_kernel_ptr;
    
    const TransformParams &param;
     
  public:

    Mat(const KernelArgs &args, const TransformParams &param) : dslash_kernel_ptr(new Kernel<KernelArgs>(args)), param(param) {}
   
    void operator()(auto &out, auto &in, auto&& transformer){
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
    
    void operator()(auto &out, auto &in){
      assert(in.GetFieldOrder() == FieldOrder::LexFieldOrder);
      // Take into account only internal points:
      const auto [Nx, Ny] = in.GetCBDims(); //Get CB dimensions

      auto X = std::views::iota(0, Nx-1);
      auto Y = std::views::iota(0, Ny-1);

      auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest
						     //
      const auto kappa = param.kappa;

      auto transformer = [=](const auto &x, const auto &y) {return (x-kappa*y);};

      auto DslashKernel = [&dslash_kernel = *dslash_kernel_ptr, 
	                   transformer_   = transformer, 
			   out_           = out.Get(), 
			   in_            = in.Get()           ] (const auto i) { 
                             //
                             dslash_kernel.apply(transformer_, out_, in_, i); 
                           };
      //
      std::for_each(std::execution::par_unseq,
                    idx.begin(),
                    idx.end(),
                    DslashKernel);
    }    
        
};


