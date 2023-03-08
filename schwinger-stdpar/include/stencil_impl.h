#pragma once
#include <algorithm>
#include <execution>
//
#include <cartesian_product.hpp>

template<typename T>
class StencilParam{
  public:
    const T kappa;	
};

template <typename gauge_tp, typename param_tp >
class StencilArgs{
  public:
    using gauge_data_tp  = typename gauge_tp::data_tp;	  

    const gauge_tp  gauge;
    
    const param_tp param;

    StencilArgs( const gauge_tp &gauge, const param_tp &param) : gauge(gauge), param(param) {}
};

template <typename Arg>
class Stencil{
  public:

    const Arg &args;

    Stencil(const Arg &args) : args(args) {}     

    void apply(auto &out, const auto &in, const auto cartesian_coords) {

      using DataTp = typename std::remove_cvref_t<Arg>::gauge_data_tp;

      auto ix = [=](auto x){ 
        return DataTp(-x.imag(), x.real());
      };

      auto [y, x] = cartesian_coords;
 
      auto o       = out.Accessor();
      const auto i = in.Accessor();      

      const auto u = args.gauge.Accessor();

      const auto kappa = args.param.kappa;
      // Update upper spin components::                       
      //mu = 0 
      auto D = u(x,y,0)*(i(x+1,y,0) - i(x+1,y,1)) + conj(u(x-1,y,0))*(i(x-1,y,0) + i(x-1,y,1)); 

      //mu = 1
      D = D + u(x,y,1)*(i(x,y+1,0) + ix(i(x,y+1,1))) + conj(u(x,y-1,1))*(i(x,y-1,0) - ix(i(x,y-1,1)));
      //                       
      o(x,y,0) = i(x,y,0) - kappa*D;
                          
      // Update down spin components::                        
      //mu = 0 
      D = u(x,y,0)*(i(x+1,y,1) - i(x+1,y,0)) + conj(u(x-1,y,0))*(i(x-1,y,1) + i(x-1,y,0));

      //mu = 1
      D = D + u(x,y,1)*(i(x,y+1,1) - ix(i(x,y+1,0))) + conj(u(x,y-1,1))*(i(x,y-1,1) + ix(i(x,y-1,0)));
                          
      o(x,y,1) = i(x,y,1) - kappa*D;        
    }
};


template<typename Kernel, typename KernelArgs>
class Mat{
  private:
    std::unique_ptr<Kernel> stencil_kernel_ptr;	

  public:

    Mat(const KernelArgs &args) : stencil_kernel_ptr(new Kernel(args)) {}

    void operator()(auto out, auto in){
       // Take into account only internal points:
       const auto [Nx, Ny] = in.GetCBDims(); //Get CB dimensions

       auto X = std::views::iota(1, Nx-2);
       auto Y = std::views::iota(1, Ny-2);

       auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest	    
       	    
       auto StencilKernel = [&stencil_kernel = *stencil_kernel_ptr.get(), out_ = out.Get(), in_ = in.Get()] (const auto i) { stencil_kernel.apply(out_, in_, i); };    
       //
       std::for_each(std::execution::par_unseq,
                     idx.begin(),
                     idx.end(),
                     StencilKernel);       
    }
};


template<typename SpinorField, typename GaugeField, typename Param>
void dispatch_dslash_kernel(SpinorField &out, GaugeField &gauge, SpinorField &in, const Param &param) {	
  assert(in.GetFieldOrder() == FieldOrder::LexFieldOrder);
  
  // Take into account only internal points:
  // Dslash_nm = (M + 4r) \delta_nm - 0.5 * \sum_\mu  ((r - \gamma_\mu)*U_(x){\mu}*\delta_{m,n+\mu} + (r + \gamma_\mu)U^*(x-mu)_{\mu}\delta_{m,n-\mu})
  //
  // gamma_{1/2} -> sigma_{1/2}, gamma_{5} -> sigma_{3}
  //
  auto &&u_ref    = gauge.Get();
  using gauge_tp  = typename std::remove_cvref_t<decltype(u_ref)>;  

  std::unique_ptr<StencilArgs<gauge_tp, Param>> stencil_args_ptr(new StencilArgs{u_ref, param});

  auto &stencil_args = *stencil_args_ptr;

  // Create dslash matrix
  auto mat = Mat<Stencil<decltype(stencil_args)>, decltype(stencil_args)>{stencil_args};
  // Apply dslash:
  mat(out, in);
}



