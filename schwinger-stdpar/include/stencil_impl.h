#pragma once
#include <algorithm>
#include <execution>
//
#include <cartesian_product.hpp>

#if 0
template<typename SpinorField, typename GaugeField, typename Arg>
void dispatch_dslash_kernel(SpinorField &out, GaugeField &gauge, SpinorField &in, const Arg &dslashArgs) {
  assert(in.GetFieldOrder() != FieldOrder::LexFieldOrder);
  
  // Take into account only internal points:
  using spinor_data_tp = typename SpinorField::data_tp;
  using gauge_data_tp  = typename GaugeField::data_tp;  

  const auto [Nx, Ny] = in.GetCBDims(); //Get CB dimensions

  auto X = std::views::iota(1, Nx-2);
  auto Y = std::views::iota(1, Ny-2);

  auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest
  // Dslash_nm = (M + 4r) \delta_nm - 0.5 * \sum_\mu  ((r - \gamma_\mu)*U_(x){\mu}*\delta_{m,n+\mu} + (r + \gamma_\mu)U^*(x-mu)_{\mu}\delta_{m,n-\mu})
  //
 // gamma_{1/2} -> sigma_{1/2}, gamma_{5} -> sigma_{3}
  const auto kappa = dslashArgs.kappa;
#if 0  						 
  // Create the kernel:
  auto DslashKernel = [=, o      = out.Accessor(), 
                          u      = gauge.Accessor(), 
                          i      = in.Accessor()](auto cartesian_coords) {
                          
                          auto ix = [=](auto x){ 
                               return spinor_data_tp(-x.imag(), x.real());
                          };
                          auto [y, x] = cartesian_coords;
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
                          
                        };
#else
  // Create the kernel:
  auto DslashKernel = [=, o_      = out, 
                          u_      = gauge, 
                          i_      = in](auto cartesian_coords) {

                          auto o      = o_.Accessor(); 
                          auto u      = u_.Accessor(); 
                          auto i      = i_.Accessor();                          
                          
                          auto ix = [=](auto x){ 
                               return spinor_data_tp(-x.imag(), x.real());
                          };
                          auto [y, x] = cartesian_coords;
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
                          
                        };
#endif                        
  std::for_each(std::execution::par_unseq,
    		idx.begin(),
    		idx.end(), 
                DslashKernel);  
}
#else

template<typename SpinorField, typename GaugeField, typename Arg>
void dispatch_dslash_kernel(SpinorField *out, GaugeField *gauge, SpinorField *in, const Arg &dslashArgs) {
  assert(in.GetFieldOrder() != FieldOrder::LexFieldOrder);
  
  // Take into account only internal points:
  using spinor_data_tp = typename SpinorField::data_tp;
  using gauge_data_tp  = typename GaugeField::data_tp;  

  const auto [Nx, Ny] = in->GetCBDims(); //Get CB dimensions

  auto X = std::views::iota(1, Nx-2);
  auto Y = std::views::iota(1, Ny-2);

  auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest
						 
  // Dslash = (I-\kappa \Sum_{mu} [ (I - \sigma_{\mu})U_(x){\mu} \delta_{x+\mu} + (I + \sigma_{\mu})U^*(x-mu)_{\mu}) \delta_{x-mu}]

  const auto kappa = dslashArgs.kappa;
#if 0  						 
  // Create the kernel:
  auto DslashKernel = [=, o      = out.Accessor(), 
                          u      = gauge.Accessor(), 
                          i      = in.Accessor()](auto cartesian_coords) {
                          
                          auto ix = [=](auto x){ 
                               return spinor_data_tp(-x.imag(), x.real());
                          };
                          auto [y, x] = cartesian_coords;
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
                          
                        };
#else
  // Create the kernel:
  auto DslashKernel = [=, o_      = out, 
                          u_      = gauge, 
                          i_      = in](auto cartesian_coords) {

                          auto o      = o_->Accessor(); 
                          auto u      = u_->Accessor(); 
                          auto i      = i_->Accessor();                          
                          
                          auto ix = [=](auto x){ 
                               return spinor_data_tp(-x.imag(), x.real());
                          };
                          auto [y, x] = cartesian_coords;
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
                          
                        };
#endif                        
  std::for_each(std::execution::par_unseq,
    		idx.begin(),
    		idx.end(), 
                DslashKernel);  
}


#endif


