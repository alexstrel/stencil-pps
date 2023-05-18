#pragma once
#include <optional>
#include <typeinfo>
//
#include <kernels/dslash.h>
#include <core/cartesian_product.hpp>

#include <typeinfo>

template<typename KernelArgs, template <typename Args> class Kernel, typename TransformParams>
class Mat{
  private:
    std::unique_ptr<Kernel<KernelArgs>> dslash_kernel_ptr;
    
    const TransformParams &param;
     
  public:

    Mat(const KernelArgs &args, const TransformParams &param) : dslash_kernel_ptr(new Kernel<KernelArgs>(args)), param(param) {}

    void operator()(SpinorFieldTp auto &out, SpinorFieldTp auto &in, auto&& transformer, const FieldOrder order = FieldOrder::EOFieldOrder){
      if ( order == FieldOrder::EOFieldOrder ) {       
        // Extract dims:
        const auto [Nx, Ny] = in.GetCBDims(); //Get CB dimensions
      
        auto X = std::views::iota(0, Nx);
        auto Y = std::views::iota(0, Ny);

        auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest
      
        auto [out_e, out_o] = out.EODecompose();
        auto [in_e,  in_o ] = in.EODecompose();       

        auto DslashEOKernel = [=, &dslash_kernel   = *dslash_kernel_ptr] (const auto coords) { 
                                //
                                dslash_kernel.template apply(transformer, out_e, in_o, coords, FieldParity::EvenFieldParity); 
                              };
        auto DslashOEKernel = [=, &dslash_kernel   = *dslash_kernel_ptr] (const auto coords) {
                                //
                                dslash_kernel.template apply(transformer, out_o, in_e, coords, FieldParity::OddFieldParity);
                              };      
        //
        std::for_each(std::execution::par_unseq,
                      idx.begin(),
                      idx.end(),
                      DslashEOKernel);

        std::for_each(std::execution::par_unseq,
                      idx.begin(),
                      idx.end(),
                      DslashOEKernel);
      } else {
        // Extract dims:      
        const auto [Nx, Ny] = in.GetDims();// in.GetCBDims(); //Get CB dimensions?
      
        auto X = std::views::iota(0, Nx);
        auto Y = std::views::iota(0, Ny);

        auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest

        auto &&out_ = out.View();
        auto &&in_  = in.View();       
        
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
    }

    void operator()(BlockSpinorFieldTp auto &out_block_spinor, BlockSpinorFieldTp auto &in_block_spinor){    
      assert(in_block_spinor.GetFieldOrder() == FieldOrder::LexFieldOrder);
      // Take into account only internal points:
      const auto [Nx, Ny] = in_block_spinor.GetDims(); //Get CB dimensions

      auto X = std::views::iota(0, Nx);
      auto Y = std::views::iota(0, Ny);

      auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest
						     //
      const auto mass = param.M;
      const auto r    = param.r;      

      const auto scale1 = mass + static_cast<decltype(r)>(2.0)*r;
      const auto scale2 = static_cast<decltype(r)>(0.5);      

      auto transformer = [=](const auto &x, const auto &y) {return (scale1*x-scale2*y);};

      //First, we need to convert to views all components in the block
      auto &&out = out_block_spinor.Convert();
      auto &&in  = in_block_spinor.Convert();       

      using spinor_ref_t =  typename std::remove_cvref_t<decltype(out.BlockView())>;

      auto &&out_ = out.BlockView();
      auto &&in_  = in.BlockView();

      auto DslashKernel = [=, &dslash_kernel = *dslash_kernel_ptr] (const auto i) mutable { 
                             //
                             dslash_kernel.template apply<spinor_ref_t>(transformer, out_, in_, i); 
                           };
      //
      std::for_each(std::execution::par_unseq,
                    idx.begin(),
                    idx.end(),
                    DslashKernel);      
    } 
};


