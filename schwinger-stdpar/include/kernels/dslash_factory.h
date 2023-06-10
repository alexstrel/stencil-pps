#pragma once
#include <optional>
#include <typeinfo>
//
#include <kernels/dslash.h>
#include <core/cartesian_product.hpp>

#include <typeinfo>

template<typename KernelArgs, template <typename Args> class Kernel>
class Mat{
  private:
    std::unique_ptr<Kernel<KernelArgs>> dslash_kernel_ptr;
     
  public:

    Mat(const KernelArgs &args) : dslash_kernel_ptr(new Kernel<KernelArgs>(args)) {}
    
    inline void launch_dslash(GenericSpinorFieldViewTp auto &out_view, GenericSpinorFieldViewTp auto &in_view, GenericSpinorFieldViewTp auto &accum_view, auto&& transformer, const FieldParity parity, const auto ids) {
      
      auto DslashKernel = [=, &dslash_kernel   = *dslash_kernel_ptr] (const auto coords) { 
                            //
                            dslash_kernel.template apply(transformer, out_view, in_view, accum_view, coords, parity); 
                          };
      //
      std::for_each(std::execution::par_unseq,
                    ids.begin(),
                    ids.end(),
                    DslashKernel);    
    }

    void operator()(GenericSpinorFieldTp auto &out, GenericSpinorFieldTp auto &in, GenericSpinorFieldTp auto &accum, auto&& transformer, const FieldParity parity){
      
      if ( in.GetFieldOrder() != FieldOrder::EOFieldOrder and in.GetFieldSubset() != FieldSiteSubset::ParitySiteSubset ) { 
        std::cerr << "Only parity field is allowed." << std::endl; 
        std::quick_exit( EXIT_FAILURE );  
      }    

      using spinor_tp    = typename std::remove_cvref_t<decltype(in)>;
      using container_tp = spinor_tp::container_tp;
      
      //Setup exe domain
      const auto [Nx, Ny] = out.GetCBDims(); //Get CB dimensions
      
      auto X = std::views::iota(0, Nx);
      auto Y = std::views::iota(0, Ny);

      auto ids = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest
      
      if constexpr (is_allocator_aware_type<container_tp>) {
        auto&& out_view    = out.View();
        auto&& in_view     = in.View();
        auto&& accum_view  = accum.View();         
        
        launch_dslash(out_view, in_view, accum_view, transformer, parity, ids);      
      } else {
        launch_dslash(out, in, accum, transformer, parity, ids);            
      }       
    }
    
    void operator()(GenericBlockSpinorFieldTp auto &out_block_spinor, GenericBlockSpinorFieldTp auto &in_block_spinor, GenericBlockSpinorFieldTp auto &accum_block_spinor, auto&& transformer, const FieldParity parity){ 
      //   
      assert(in_block_spinor.GetFieldOrder() == FieldOrder::EOFieldOrder and in_block_spinor.GetFieldSubset() == FieldSiteSubset::ParitySiteSubset);
      
      using block_spinor_tp        = typename std::remove_cvref_t<decltype(in_block_spinor)>;
      using component_container_tp = block_spinor_tp::container_tp;      
      
      // Take into account only internal points:
      const auto [Nx, Ny] = in_block_spinor.GetCBDims(); //Get CB dimensions

      auto X = std::views::iota(0, Nx);
      auto Y = std::views::iota(0, Ny);

      auto ids = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest
                    
     if constexpr (is_allocator_aware_type<component_container_tp>) {
        //First, we need to convert to views all components in the block
        auto &&out_block_spinor_view    = out_block_spinor.ConvertToView();
        auto &&in_block_spinor_view     = in_block_spinor.ConvertToView();       
        auto &&accum_block_spinor_view  = accum_block_spinor.ConvertToView();               

        auto &&out_view    = out_block_spinor_view.BlockView();
        auto &&in_view     = in_block_spinor_view.BlockView(); 
        auto &&accum_view  = accum_block_spinor_view.BlockView();         
        
        launch_dslash(out_view, in_view, accum_view, transformer, parity, ids);      
      } else {
        auto &&out_view    = out_block_spinor.BlockView();
        auto &&in_view     = in_block_spinor.BlockView(); 
        auto &&accum_view  = accum_block_spinor.BlockView();         
      
        launch_dslash(out_view, in_view, accum_view, transformer, parity, ids);      
      }                    
    } 
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
};


