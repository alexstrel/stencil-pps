#pragma once

template<typename SpinorField, typename GaugeField>
void dispatch_dslash_kernel(SpinorField &out, const GaugeField &gauge, const SpinorField &in) {
  // Take into account only internal points:
  const auto [Nx, Ny] = in.GetDims(); //Get CB dimensions

  auto xh = std::views::iota(0, in.GetFieldSubset == FieldSubset::FullSiteSubset ? (Nx / 2 - 1) : (Nx   -1));
  auto y  = std::views::iota(0, Ny -1);

  auto idx = std::views::cartesian_product(y, xh);//Y is the slowest index, X is the fastest
						 
  auto [even_in, odd_in]   = in.EODecompose();
  auto [even_out, odd_out] = out.EODecompose();

						 
  // Create the kernel:
  auto stencil_kernel = [out = out.data(), in = in.data(), args](auto cartesian_coords) {
                          
                          auto idx = [=](auto x, auto y){ 
                               return x + y*args.nx;
                          };
                          
                          auto [y, x] = cartesian_coords;
                                                               
                          // Load hopping terms or apply boundary conditions:                         
                          const auto in_xm1 = x == 1             ? 0 : in[idx(x - 1, y)];
                          const auto in_xp1 = x == (args.nx - 2) ? 0 : in[idx(x + 1, y)];                          
                          //
                          const auto in_ym1 = y == 1             ? 1 : in[idx(x, y - 1)];
                          const auto in_yp1 = y == (args.ny - 2) ? 0 : in[idx(x, y + 1)]; 

                          auto site_stencil = (args.C(0)*in[idx(x, y)] + args.C(1) * (in_ym1 + in_yp1 + in_xm1 + in_xp1));

                          //store the result
                          out[idx(x, y)] = site_stencil;
                          // we may want to return something if we plan to run the reduction-like algorithm
                          return compute_energy ? site_stencil*args.dx*args.dx : 0.0;
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

