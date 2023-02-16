#pragma once

template<typename SpinorField, typename GaugeField>
void dispatch_dslash_kernel(SpinorField &out, const GaugeField &gauge, const SpinorField &in) {
  // Take into account only internal points:
  const auto [Nxh, Ny] = in.GetCBDims(); //Get CB dimensions

  auto xh = std::views::iota(0, Nxh-1);
  auto y  = std::views::iota(0, Ny -1);

  auto idx = std::views::cartesian_product(y, xh);//Y is the slowest index, X is the fastest
						 
  auto [even_out, odd_out] = out.EODecompose();

  const auto [even_in, odd_in]       = in.EODecompose();
  const auto [even_gauge, odd_gauge] = gauge.EODecompose();  

  // Dslash = (I-\kappa \Sum_{mu} [ (I - \sigma_{\mu})U_(x){\mu} \delta_{x+\mu} + (I + \sigma_{\mu})U^*(x-mu)_{\mu}) \delta_{x-mu}]
  //        = I_{EE}  D_{EO}
  //          D_{OE}  I_{OO}
						 
  // Create the kernel:
  auto Dslash = [&o_e      = even_out->Accessor(), 
                 &o_o      = odd_out->Accessor(), 
                 &U_e      = even_gauge->Accessor(), 
                 &U_o      = odd_gauge->Accessor(), 
                 &i_e      = even_in->Accessor(), 
                 &i_o      = odd_in->Accessor()](auto cartesian_coords) {
                          
                          auto [y_, xh_] = cartesian_coords;
                          
                          const auto other_parity_site_forward_hop = y_ & 1; //(y+z+t) & 1 => parity index
                          const auto other_parity_site_bckward_hop = other_site_forward_hop - 1;                           
                          
                          const auto xhp1  = xh_ + other_site_forward_hop;
                          const auto xhm1  = xh_ + other_site_bckward_hop; 
                          //
                          const auto yp1   = y_ + 1;
                          const auto ym1   = y_ - 1; 
                           
                          auto site_stencil = (args.C(0)*in[idx(x, y)] + args.C(1) * (in_ym1 + in_yp1 + in_xm1 + in_xp1));

                          //store the result
                          out[idx(x, y)] = site_stencil;
                          // we may want to return something if we plan to run the reduction-like algorithm
                        };
                        
  if constexpr (compute_energy){                          
    auto energy = std::transform_reduce(std::execution::par_unseq,
    					idx.begin(), 
    					idx.end(),  
                                        0., 
                                        std::plus{}, 
                                        stencil_kernel);  
                                        
    return std::optional<FloatTp>(energy);                                      
  } else {
    std::for_each(std::execution::par_unseq,
    		  idx.begin(),
    		  idx.end(), 
                  stencil_kernel);  

    return std::nullopt;
  }
}

