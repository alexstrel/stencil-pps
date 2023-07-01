#pragma once
#include <optional>
#include <typeinfo>
//
#include <kernels/dslash.h>
#include <core/cartesian_product.hpp>
//
#include <fields/field.h>

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
    
    template<bool dagger>
    inline void launch_dslash(GenericSpinorFieldViewTp auto &out_view, const GenericSpinorFieldViewTp auto &in_view, const GenericSpinorFieldViewTp auto &aux_view, auto&& post_transformer, const FieldParity parity, const auto ids) {
      
      auto DslashKernel = [=, &dslash_kernel   = *dslash_kernel_ptr] (const auto coords) { 
                            //
                            dslash_kernel.template apply<dagger>(out_view, in_view, aux_view, post_transformer, coords, parity); 
                          };
      //
      std::for_each(std::execution::par_unseq,
                    ids.begin(),
                    ids.end(),
                    DslashKernel);    
    } 
    
    template<bool dagger>
    inline void launch_dslash(GenericSpinorFieldViewTp auto &out_view, const GenericSpinorFieldViewTp auto &in_view, const FieldParity parity, const auto ids) {
      
      auto DslashKernel = [=, &dslash_kernel   = *dslash_kernel_ptr] (const auto coords) { 
                            //
                            dslash_kernel.template apply<dagger>(out_view, in_view, coords, parity); 
                          };
      //
      std::for_each(std::execution::par_unseq,
                    ids.begin(),
                    ids.end(),
                    DslashKernel);    
    }    
 
    void operator()(GenericSpinorFieldTp auto &out, const GenericSpinorFieldTp auto &in, const GenericSpinorFieldTp auto &aux, auto&& post_transformer, const FieldParity parity, const bool dagger = false){
      
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
        auto&& out_view       = out.View();
        const auto&& in_view  = in.View();
        const auto&& aux_view = aux.View();         
        
        if (dagger) {
          launch_dslash<true>(out_view, in_view, aux_view, post_transformer, parity, ids);
        } else {
          launch_dslash<false>(out_view, in_view, aux_view, post_transformer, parity, ids);
        }      
      } else {
        if (dagger) {
          launch_dslash<true>(out, in, aux, post_transformer, parity, ids);  
        } else {
          launch_dslash<false>(out, in, aux, post_transformer, parity, ids);
        }
      }       
    }
    
    void operator()(GenericBlockSpinorFieldTp auto &out_block_spinor, GenericBlockSpinorFieldTp auto &in_block_spinor, GenericBlockSpinorFieldTp auto &aux_block_spinor, auto&& post_transformer, const FieldParity parity, const bool dagger = false){ 
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
        auto &&aux_block_spinor_view  = aux_block_spinor.ConvertToView();               

        auto &&out_view    = out_block_spinor_view.BlockView();
        auto &&in_view     = in_block_spinor_view.BlockView(); 
        auto &&aux_view  = aux_block_spinor_view.BlockView();         
        
        if (dagger) {
          launch_dslash<true>(out_view, in_view, aux_view, post_transformer, parity, ids);  
        } else {
          launch_dslash<false>(out_view, in_view, aux_view, post_transformer, parity, ids);
        }  
      } else {
        auto &&out_view    = out_block_spinor.BlockView();
        auto &&in_view     = in_block_spinor.BlockView(); 
        auto &&aux_view  = aux_block_spinor.BlockView();         
      
        if (dagger) {
          launch_dslash<true>(out_view, in_view, aux_view, post_transformer, parity, ids);  
        } else {
          launch_dslash<false>(out_view, in_view, aux_view, post_transformer, parity, ids);
        }     
      }                    
    }    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
};

template<typename KernelArgs, template <typename Args> class Kernel, typename TransformParams>
class Mat : public MatTransform<KernelArgs, Kernel> {
  private:
    const TransformParams param;
    
    const FieldParity parity;
 
    bool base_dagger;    
  public:

    Mat(const KernelArgs &args, const TransformParams &param, const FieldParity parity = FieldParity::InvalidFieldParity, const bool dagger = false) : MatTransform<KernelArgs, Kernel>(args), param(param), parity(parity), base_dagger(dagger) {}
    
    inline void flip() { base_dagger = not base_dagger; }

    void operator()(GenericSpinorFieldTp auto &out, const GenericSpinorFieldTp auto &in, const GenericSpinorFieldTp auto &aux){
      // Check all arguments!
      if(out.GetFieldSubset() != FieldSiteSubset::ParitySiteSubset) { 
        std::cerr << "Error: undefined parity.. exiting\n";
        std::quick_exit( EXIT_FAILURE );
      }       
      
      const auto const1 = static_cast<MatTransform<KernelArgs, Kernel>::kernel_data_tp>(param.M + 2.0*param.r);
      const auto const2 = static_cast<MatTransform<KernelArgs, Kernel>::kernel_data_tp>(0.5);                
      
      auto transformer = [=](const auto &x, const auto &y) {return (const1*x-const2*y);};      
      //
      MatTransform<KernelArgs, Kernel>::operator()(out, in,  aux, transformer, parity, base_dagger);
    }
    
    void operator()(GenericSpinorFieldTp auto &out, GenericSpinorFieldTp auto &in){//FIXME: in argument must be constant
      // Check all arguments!
      if( out.GetFieldSubset() != FieldSiteSubset::FullSiteSubset ) { 
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
      MatTransform<KernelArgs, Kernel>::operator()(even_out, odd_in,  even_in, transformer, FieldParity::EvenFieldParity, base_dagger);
      MatTransform<KernelArgs, Kernel>::operator()(odd_out,  even_in, odd_in,  transformer, FieldParity::OddFieldParity, base_dagger);       
    }    
    
    void operator()(GenericBlockSpinorFieldTp auto &out_block_spinor, const GenericBlockSpinorFieldTp auto &in_block_spinor, const GenericBlockSpinorFieldTp auto &aux_block_spinor){ 
      //
      // Check all arguments!
      if( out_block_spinor.GetFieldSubset() != FieldSiteSubset::ParitySiteSubset ) { 
        std::cerr << "Error: undefined parity.. exiting\n";
        std::quick_exit( EXIT_FAILURE );
      }       
      
      const auto const1 = static_cast<MatTransform<KernelArgs, Kernel>::kernel_data_tp>(param.M + 2.0*param.r);
      const auto const2 = static_cast<MatTransform<KernelArgs, Kernel>::kernel_data_tp>(0.5);                
      
      auto transformer = [=](const auto &x, const auto &y) {return (const1*x-const2*y);};      
      //
      MatTransform<KernelArgs, Kernel>::operator()(out_block_spinor, in_block_spinor, aux_block_spinor, transformer, parity, base_dagger);         
    }   
    
 void operator()(GenericBlockSpinorFieldTp auto &out_block_spinor, GenericBlockSpinorFieldTp auto &in_block_spinor){ 
      //
      // Check all arguments!
      if(out_block_spinor.GetFieldSubset() != FieldSiteSubset::FullSiteSubset) { 
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
      MatTransform<KernelArgs, Kernel>::operator()(even_out_block, odd_in_block,  even_in_block, transformer, FieldParity::EvenFieldParity, base_dagger);
      MatTransform<KernelArgs, Kernel>::operator()(odd_out_block,  even_in_block, odd_in_block,  transformer, FieldParity::OddFieldParity, base_dagger);       
    }    
     
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
};

template<typename KernelArgs, template <typename Args> class Kernel, typename TransformParams, typename Spinor, bool do_normal = false>
class PreconMat : public MatTransform<KernelArgs, Kernel> {
  private:
    const TransformParams param;
    
    Spinor tmp;    
    Spinor tmp2;
    
    const FieldParity parity;
 
    bool base_dagger; 
       
  public:

    using ArgTp        = decltype(std::declval<Spinor>().ExportArg());
    using container_tp = Spinor::container_tp;

    PreconMat(const KernelArgs &args, const TransformParams &param, const Spinor &spinor,  const FieldParity parity = FieldParity::InvalidFieldParity, const bool dagger = false) : 
                                            MatTransform<KernelArgs, Kernel>(args), 
                                            param(param), 
                                            tmp{create_field<container_tp, ArgTp>(spinor.ExportArg())},
                                            tmp2{create_field<container_tp, ArgTp>(spinor.ExportArg())},
                                            parity(parity), 
                                            base_dagger(dagger) { 
                                              assert(tmp.GetFieldSubset() == FieldSiteSubset::ParitySiteSubset);
                                            }
   
    
    inline void flip() { base_dagger = not base_dagger; }

    void operator()(GenericSpinorFieldTp auto &out, GenericSpinorFieldTp auto &in){//FIXME: in argument must be constant
      // Check all arguments!
      if(out.GetFieldSubset() != FieldSiteSubset::ParitySiteSubset) { 
        std::cerr << "This operation is supported for parity fields only...\n";
        std::quick_exit( EXIT_FAILURE );
      }            
      const auto c = static_cast<MatTransform<KernelArgs, Kernel>::kernel_data_tp>( 1.0 / (2.0 * ( param.M + 2.0*param.r)));
      
      auto transformer = [=](const auto &x, const auto &y) {return (x - c*y);};      
      //
      auto other_parity = parity == EvenFieldParity ? FieldParity::OddFieldParity : FieldParity::EvenFieldParity;
      //
      MatTransform<KernelArgs, Kernel>::operator()(tmp,  in, parity, base_dagger);
      MatTransform<KernelArgs, Kernel>::operator()(out,  tmp, in,  transformer, other_parity, base_dagger);       

      if constexpr (do_normal) {
        MatTransform<KernelArgs, Kernel>::operator()(tmp,  in, parity, base_dagger);
        MatTransform<KernelArgs, Kernel>::operator()(tmp2,  tmp, in,  transformer, other_parity, base_dagger);

        flip();
        MatTransform<KernelArgs, Kernel>::operator()(tmp,  tmp2, parity, base_dagger);
        MatTransform<KernelArgs, Kernel>::operator()(out,  tmp, tmp2,  transformer, other_parity, base_dagger);
        flip();
      } else {
        MatTransform<KernelArgs, Kernel>::operator()(tmp,  in, parity, base_dagger);
        MatTransform<KernelArgs, Kernel>::operator()(out,  tmp, in,  transformer, other_parity, base_dagger);
      }
    }    
    
};

