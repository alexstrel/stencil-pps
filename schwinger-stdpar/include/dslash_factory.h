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
      const auto [Nx, Ny] = in.GetDims();// in.GetCBDims(); //Get CB dimensions?
      
      auto X = std::views::iota(0, Nx);
      auto Y = std::views::iota(0, Ny);

      auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest
      
      auto &&out_ = out.Reference();
      auto &&in_  = in.Reference();       

      auto DslashKernel = [=, &dslash_kernel   = *dslash_kernel_ptr] (const auto coords) { 
                                //
                                dslash_kernel.template apply(transformer, out_, in_, coords); 
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
      const auto [Nx, Ny] = in.GetDims(); //Get CB dimensions

      auto X = std::views::iota(0, Nx);
      auto Y = std::views::iota(0, Ny);

      auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest
						     //
      const auto mass = param.M;
      const auto r    = param.r;      

      auto transformer = [=](const auto &x, const auto &y) {return ((mass+2.0*r)*x-0.5*y);};

      auto &&out_ = out.Reference();
      auto &&in_  = in.Reference();       

      auto DslashKernel = [=, &dslash_kernel = *dslash_kernel_ptr] (const auto i) { 
                             //
                             dslash_kernel.template apply(transformer, out_, in_, i); 
                           };
      //
      std::for_each(std::execution::par_unseq,
                    idx.begin(),
                    idx.end(),
                    DslashKernel);
    }    
};


