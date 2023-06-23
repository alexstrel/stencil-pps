#pragma once
#include <optional>
#include <typeinfo>
//
#include <kernels/dslash.h>
#include <core/cartesian_product.hpp>

#include <typeinfo>

//MatTransform
template<typename KernelArgs, template <typename Args> class Kernel>
class MatTransform{
  private:
    std::unique_ptr<Kernel<KernelArgs>> dslash_kernel_ptr;
     
  public:
    using kernel_data_tp = typename std::remove_cvref_t<KernelArgs>::gauge_data_tp;

    MatTransform(const KernelArgs &args) : dslash_kernel_ptr(new Kernel<KernelArgs>(args)) {}
    
    KernelArgs& ExportKernelArgs() const { return dslash_kernel_ptr->args; }
    
    inline void launch_dslash(GenericSpinorFieldViewTp auto &out_view, const GenericSpinorFieldViewTp auto &in_view, const GenericSpinorFieldViewTp auto &accum_view, auto&& post_transformer, const FieldParity parity, const auto ids) {
      
      auto DslashKernel = [=, &dslash_kernel   = *dslash_kernel_ptr] (const auto coords) { 
                            //
                            dslash_kernel.template apply(out_view, in_view, accum_view, post_transformer, coords, parity); 
                          };
      //
      std::for_each(std::execution::par_unseq,
                    ids.begin(),
                    ids.end(),
                    DslashKernel);    
    }

    void operator()(GenericSpinorFieldTp auto &out, const GenericSpinorFieldTp auto &in, const GenericSpinorFieldTp auto &accum, auto&& post_transformer, const FieldParity parity){
      
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
        const auto&& in_view     = in.View();
        const auto&& accum_view  = accum.View();         
        
        launch_dslash(out_view, in_view, accum_view, post_transformer, parity, ids);      
      } else {
        launch_dslash(out, in, accum, post_transformer, parity, ids);            
      }       
    }
    
    void operator()(GenericBlockSpinorFieldTp auto &out_block_spinor, GenericBlockSpinorFieldTp auto &in_block_spinor, GenericBlockSpinorFieldTp auto &accum_block_spinor, auto&& post_transformer, const FieldParity parity){ 
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
        
        launch_dslash(out_view, in_view, accum_view, post_transformer, parity, ids);      
      } else {
        auto &&out_view    = out_block_spinor.BlockView();
        auto &&in_view     = in_block_spinor.BlockView(); 
        auto &&accum_view  = accum_block_spinor.BlockView();         
      
        launch_dslash(out_view, in_view, accum_view, post_transformer, parity, ids);      
      }                    
    }    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
};

template<typename KernelArgs, template <typename Args> class Kernel, typename TransformParams>
class Mat : public MatTransform<KernelArgs, Kernel> {
  private:
    const TransformParams param;
    
    const FieldParity parity;
     
  public:

    Mat(const KernelArgs &args, const TransformParams &param, const FieldParity parity = FieldParity::InvalidFieldParity) : MatTransform<KernelArgs, Kernel>(args), param(param), parity(parity) {}

    void operator()(GenericSpinorFieldTp auto &out, const GenericSpinorFieldTp auto &in, const GenericSpinorFieldTp auto &accum){
      // Check all arguments!
      if(parity == FieldParity::InvalidFieldParity) { 
        std::cerr << "Error: undefined parity.. exiting\n";
        std::quick_exit( EXIT_FAILURE );
      }       
      
      const auto const1 = static_cast<MatTransform<KernelArgs, Kernel>::kernel_data_tp>(param.M + 2.0*param.r);
      const auto const2 = static_cast<MatTransform<KernelArgs, Kernel>::kernel_data_tp>(0.5);                
      
      auto transformer = [=](const auto &x, const auto &y) {return (const1*x-const2*y);};      
      //
      MatTransform<KernelArgs, Kernel>::operator()(out, in,  accum, transformer, parity);
    }
    
    void operator()(GenericSpinorFieldTp auto &out, GenericSpinorFieldTp auto &in){//FIXME: in argument must be constant
      // Check all arguments!
      if(parity != FieldParity::InvalidFieldParity) { 
        std::cerr << "This operation is supported for full fields only...\n";
        std::quick_exit( EXIT_FAILURE );
      }            
      const auto const1 = static_cast<MatTransform<KernelArgs, Kernel>::kernel_data_tp>(param.M + 2.0*param.r);
      const auto const2 = static_cast<MatTransform<KernelArgs, Kernel>::kernel_data_tp>(0.5);                
      
      auto transformer = [=](const auto &x, const auto &y) {return (const1*x-const2*y);};      
      //
      auto [even_in,   odd_in] = in.EODecompose();
      auto [even_out, odd_out] = out.EODecompose();      
      //
      MatTransform<KernelArgs, Kernel>::operator()(even_out, odd_in,  even_in, transformer, FieldParity::EvenFieldParity);
      MatTransform<KernelArgs, Kernel>::operator()(odd_out,  even_in, odd_in,  transformer, FieldParity::OddFieldParity);       
    }    
    
    void operator()(GenericBlockSpinorFieldTp auto &out_block_spinor, const GenericBlockSpinorFieldTp auto &in_block_spinor, const GenericBlockSpinorFieldTp auto &accum_block_spinor){ 
      //
      // Check all arguments!
      if(parity == FieldParity::InvalidFieldParity) { 
        std::cerr << "Error: undefined parity.. exiting\n";
        std::quick_exit( EXIT_FAILURE );
      }       
      
      const auto const1 = static_cast<MatTransform<KernelArgs, Kernel>::kernel_data_tp>(param.M + 2.0*param.r);
      const auto const2 = static_cast<MatTransform<KernelArgs, Kernel>::kernel_data_tp>(0.5);                
      
      auto transformer = [=](const auto &x, const auto &y) {return (const1*x-const2*y);};      
      //
      MatTransform<KernelArgs, Kernel>::operator()(out_block_spinor, in_block_spinor, accum_block_spinor, transformer, parity);         
    }   
    
 void operator()(GenericBlockSpinorFieldTp auto &out_block_spinor, GenericBlockSpinorFieldTp auto &in_block_spinor){ 
      //
      // Check all arguments!
      if(parity != FieldParity::InvalidFieldParity) { 
        std::cerr << "This operation is supported for full fields only...\n";
        std::quick_exit( EXIT_FAILURE );
      }      
      
      const auto const1 = static_cast<MatTransform<KernelArgs, Kernel>::kernel_data_tp>(param.M + 2.0*param.r);
      const auto const2 = static_cast<MatTransform<KernelArgs, Kernel>::kernel_data_tp>(0.5);                
      
      auto transformer = [=](const auto &x, const auto &y) {return (const1*x-const2*y);};      
      //
      auto [even_in_block,   odd_in_block] = in_block_spinor.EODecompose();
      auto [even_out_block, odd_out_block] = out_block_spinor.EODecompose();      
      //
      MatTransform<KernelArgs, Kernel>::operator()(even_out_block, odd_in_block,  even_in_block, transformer, FieldParity::EvenFieldParity);
      MatTransform<KernelArgs, Kernel>::operator()(odd_out_block,  even_in_block, odd_in_block,  transformer, FieldParity::OddFieldParity);       
    }    
     
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
};


