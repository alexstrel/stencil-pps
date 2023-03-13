#pragma once
#include <dslash.h>

template<typename Kernel, typename KernelArgs, typename TransformParams>
class Mat{
  private:
    std::unique_ptr<Kernel> dslash_kernel_ptr;

    const TransformParams &param;

  public:

    Mat(const KernelArgs &args, const TransformParams &param) : dslash_kernel_ptr(new Kernel(args)), 
	                                                        param(param) {}

    void operator()(auto &out, auto &in){
      assert(in.GetFieldOrder() == FieldOrder::LexFieldOrder);
      //
      constexpr bool do_pre_transform = false;
      // Extract dims:
      const auto [Nx, Ny] = in.GetCBDims(); //Get CB dimensions

      auto X = std::views::iota(0, Nx-1);
      auto Y = std::views::iota(0, Ny-1);

      auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest
						     //
      const auto kappa = param.kappa;

      auto transformer = [=](const auto &x, const auto &y) {return (x-kappa*y);};

      auto DslashKernel = [&dslash_kernel   = *dslash_kernel_ptr, 
	                   post_transformer = transformer, 
			   out_             = out.Get(), 
			   in_              = in.Get()           ] (const auto coords) { 
                             //
                             dslash_kernel.template apply<do_pre_transform>(0, 
					                                    post_transformer, 
									    out_, 
									    in_, 
									    coords); 
                           };
      //
      std::for_each(std::execution::par_unseq,
                    idx.begin(),
                    idx.end(),
                    DslashKernel);
    }
};


