#pragma once
#include <optional>
#include <typeinfo>
//
#include <kernels/dslash.h>
#include <core/cartesian_product.hpp>

template<typename KernelArgs, template <typename Args> class Kernel, typename TransformParams>
class Mat{
  private:
    std::unique_ptr<Kernel<KernelArgs>> dslash_kernel_ptr;
    
    const TransformParams &param;
     
  public:

    Mat(const KernelArgs &args, const TransformParams &param) : dslash_kernel_ptr(new Kernel<KernelArgs>(args)), param(param) {}

    template<SpinorFieldTp spinor_tp>
    void operator()(spinor_tp &out, spinor_tp &in, auto&& transformer, const FieldOrder order = FieldOrder::EOFieldOrder){	         
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

      const auto scale1 = mass + static_cast<decltype(r)>(2.0)*r;
      const auto scale2 = static_cast<decltype(r)>(0.5);      

      auto transformer = [=](const auto &x, const auto &y) {return (scale1*x-scale2*y);};

      auto &&out_ = out.View();
      auto &&in_  = in.View();       

      using spinor_ref_t =  typename std::remove_cvref_t<decltype(out_)>;

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
/*
    template<BlockSpinorFieldTp block_spinor_field>
    void operator()(block_spinor_field &out, block_spinor_field &in){
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

      auto &&out_ = out.View();
      auto &&in_  = in.View();

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
*/    
};


