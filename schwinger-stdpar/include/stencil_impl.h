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

template <typename spinor_tp, typename gauge_tp, typename param_tp >
class StencilArgs{
  public:
    using spinor_data_tp = typename spinor_tp::data_tp;
    using gauge_data_tp  = typename gauge_tp::data_tp;	  

    spinor_tp out;
    const spinor_tp in;
    const gauge_tp  gauge;
    
    const param_tp param;

    StencilArgs( spinor_tp &out, const gauge_tp &gauge, const spinor_tp &in, const param_tp &param) : out(out), in(in), gauge(gauge), param(param) {}
};

template <typename Arg>
class Stencil{
  public:

    const Arg &args;

    Stencil(const Arg &args) : args(args) {}     

    void operator()(auto cartesian_coords) {

      using DataTp = typename std::remove_cvref_t<Arg>::spinor_data_tp;

      auto ix = [=](auto x){ 
        return DataTp(-x.imag(), x.real());
      };

      auto [y, x] = cartesian_coords;
 
      auto o      = args.out.Accessor();
      auto u      = args.gauge.Accessor();
      auto i      = args.in.Accessor();

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


template<typename SpinorField, typename GaugeField, typename Param>
void dispatch_dslash_kernel(SpinorField &out, GaugeField &gauge, SpinorField &in, const Param &param) {	
  assert(in.GetFieldOrder() == FieldOrder::LexFieldOrder);
  
  // Take into account only internal points:
  const auto [Nx, Ny] = in.GetCBDims(); //Get CB dimensions

  auto X = std::views::iota(1, Nx-2);
  auto Y = std::views::iota(1, Ny-2);

  auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest
  // Dslash_nm = (M + 4r) \delta_nm - 0.5 * \sum_\mu  ((r - \gamma_\mu)*U_(x){\mu}*\delta_{m,n+\mu} + (r + \gamma_\mu)U^*(x-mu)_{\mu}\delta_{m,n-\mu})
  //
  // gamma_{1/2} -> sigma_{1/2}, gamma_{5} -> sigma_{3}
  //
  auto &&in_ref  = in.Get();
  auto &&out_ref = out.Get();
  auto &&u_ref   = gauge.Get();

  using spinor_tp = typename std::remove_cvref_t<decltype(in_ref)>;  
  using gauge_tp  = typename std::remove_cvref_t<decltype(u_ref)>;  

  std::unique_ptr<StencilArgs<spinor_tp, gauge_tp, Param>> stencil_args_ptr(new StencilArgs{out_ref, u_ref, in_ref, param});

  auto &stencil_args = *stencil_args_ptr;

  auto dslash_kernel_ptr = std::make_shared<Stencil<decltype(stencil_args)>>(stencil_args);

  auto DslashKernel = [&dslash_kernel = *dslash_kernel_ptr.get()] (const auto i) { dslash_kernel(i); }; 

  std::for_each(std::execution::par_unseq,
    		idx.begin(),
    		idx.end(), 
                DslashKernel);  
}



