#pragma once

template<typename SpinorField, typename GaugeField>
void dispatch_dslash_kernel(SpinorField &out, const GaugeField &gauge, const SpinorField &in) {
  // Take into account only internal points:
  using spinor_data_tp = typename SpinorField::data_tp;
  using gauge_data_tp  = typename GaugeField::data_tp;  

  const auto [Nxh, Ny] = in.GetCBDims(); //Get CB dimensions

  auto Xh = std::views::iota(1, Nxh-2);
  auto Y  = std::views::iota(1, Ny -2);

  auto idx = std::views::cartesian_product(Y, Xh);//Y is the slowest index, X is the fastest
						 
  auto [even_out, odd_out] = out.EODecompose();

  const auto [even_in, odd_in]       = in.EODecompose();
  const auto [even_gauge, odd_gauge] = gauge.EODecompose();  

  // Dslash = (I-\kappa \Sum_{mu} [ (I - \sigma_{\mu})U_(x){\mu} \delta_{x+\mu} + (I + \sigma_{\mu})U^*(x-mu)_{\mu}) \delta_{x-mu}]
  //        = I_{EE}  D_{EO}
  //          D_{OE}  I_{OO}
						 
  // Create the kernel:
  auto Dslash = [&oe      = even_out->Accessor(), 
                 &oo      = odd_out->Accessor(), 
                 &ue      = even_gauge->Accessor(), 
                 &uo      = odd_gauge->Accessor(), 
                 &ie      = even_in->Accessor(), 
                 &io      = odd_in->Accessor()](auto cartesian_coords) {
                          
                          auto [y, xh] = cartesian_coords;
			
                          // Update even sites:
			  //
                          auto other_parity_bit = y & 1; //(y+z+t) & 1 => parity index
			  //
                          auto forward_hop = other_parity_bit; 							  
                          auto bckward_hop = 1 - other_parity_bit;                           
                          
                          auto xhp1 = xh + forward_hop;
                          auto xhm1 = xh - bckward_hop; 
                          //
                          const auto yp1 = y + 1;
                          const auto ym1 = y - 1; 

                          // Update upper spin components::			  
                          //mu = 0 
			  auto Deo = ue(xh,y,0)*(io(xhp1,y,0) - io(xhp1,y,1)) + conj(uo(xhm1,y,0))*(io(xhm1,y,0) + io(xhm1,y,1)); 

                          //mu = 1
			  auto io_yp1   = io(xh,yp1,1);
			  auto ixio_yp1 = spinor_data_tp(-io_yp1.imag(), io_yp1.real());

			  auto io_ym1   = io(xh,ym1,1);
			  auto ixio_ym1 = spinor_data_tp(-io_ym1.imag(), io_ym1.real()); 

                          Deo = Deo + ue(xh,y,1)*(io(xh,yp1,0) + ixio_yp1) + conj(uo(xh,ym1,1))*(io(xh,ym1,0) - ixio_ym1);
                          //                       
                          oe(xh, y, 0) = ie(xh, y, 0) - kappa*Deo;
			  
                          // Update down spin components::			  
                          //mu = 0 
                          Deo = ue(xh,y,0)*(io(xhp1,y,1) - io(xhp1,y,0)) + conj(uo(xhm1,y,0))*(io(xhm1,y,1) + io(xhm1,y,0));

                          //mu = 1
                          io_yp1   = io(xh,yp1,0);
                          ixio_yp1 = spinor_data_tp(-io_yp1.imag(), io_yp1.real());

                          io_ym1   = io(xh,ym1,0);
                          ixio_ym1 = spinor_data_tp(-io_ym1.imag(), io_ym1.real());

                          Deo = Deo + ue(xh,y,1)*(io(xh,yp1,1) - ixio_yp1) + conj(uo(xh,ym1,1))*(io(xh,ym1,1) + ixio_ym1);
			  
                          oe(xh, y, 1) = ie(xh, y, 1) - kappa*Doe;
			  
                          // Update odd sites:
			  auto other_parity_bit = 1 - other_parity_bit;
			  //
                          auto forward_hop = other_parity_bit;                                      
                          auto bckward_hop = 1 - other_parity_bit;

                          auto xhp1  = xh + forward_hop;       
                          auto xhm1  = xh - bckward_hop;

                          
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

