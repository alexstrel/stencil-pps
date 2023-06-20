#pragma once
#include <algorithm>
#include <execution>
//
#include <core/cartesian_product.hpp>
#include <kernels/dslash_helpers.h>
//
#include <core/color_spinor.h>

template<typename T>
class DslashParam{
  public:
    const T M;
    const T r;    
};


//template <GaugeFieldViewTp gauge_tp, int nSpin_ = 2, int bSize_  = 1>
template <GaugeFieldViewTp gauge_tp, int nSpin_ = 2>
class DslashArgs {
  public:
    using gauge_data_tp  = typename gauge_tp::data_tp;	  

    static constexpr std::size_t nDir   = gauge_tp::nDir;
    static constexpr std::size_t nColor = gauge_tp::nColor;
    static constexpr std::size_t nSpin  = nSpin_;
//    static constexpr std::size_t bSize  = bSize_;    
    static constexpr std::size_t bSize  = 1;    
    
    using SpinorTp = impl::ColorSpinor<gauge_data_tp, nColor, nSpin, bSize>;
    using LinkTp   = impl::ColorMatrix<gauge_data_tp, nColor, bSize>;       

    const gauge_tp  gauge;
    
    DslashArgs( const gauge_tp &gauge) : gauge(gauge) {}        
    
    
    template<std::size_t... Idxs>
    static inline decltype(auto) load_spinor(std::index_sequence<Idxs...>, const auto& field_accessor, const std::array<int, nDir>& x){
      std::array<gauge_data_tp, nColor*nSpin*bSize> spinor{field_accessor(x[Idxs]..., 0), field_accessor(x[Idxs]..., 1)};    
      return SpinorTp(spinor);//2-component spinor
    }    
        
    template<std::size_t... Idxs>
    static inline decltype(auto) load_parity_link(std::index_sequence<Idxs...>, const auto& field_accessor, const std::array<int, nDir>& x, const int &d, const int &parity){
      std::array<gauge_data_tp, nColor*nColor*bSize> link{field_accessor(x[Idxs]..., d, parity)};
      return LinkTp(link);
    }        
};

template <typename Arg>
class Dslash{
  public:
    using ArgTp  = typename std::remove_cvref_t<Arg>;

    static constexpr std::size_t nSpin = ArgTp::nSpin;
    static constexpr std::size_t nDir  = ArgTp::nDir;    

    using Link   = ArgTp::LinkTp; 
    using Spinor = ArgTp::SpinorTp;
    
    using Indices = std::make_index_sequence<nDir>;    

    const Arg &args;

    Dslash(const Arg &args) : args(args) {}     


    inline decltype(auto) compute_parity_site_stencil(const auto &in_accessor, const auto &U_accessor, const FieldParity parity, const std::array<int, nDir> site_coords){
      //Define accessor wrappers:

      auto in = [&in_=in_accessor, this](const std::array<int, nDir> &x){ 
        return ArgTp::load_spinor(Indices{}, in_, x);
      };

      auto U  = [&U_=U_accessor, this](const std::array<int, nDir> &x, const int &d, const int &p){ 
        return ArgTp::load_parity_link(Indices{}, U_, x, d, p);
      };         
    
      auto is_local_boundary = [](const auto d, const auto coord, const auto bndry, const auto parity_bit){ 
	      return ((coord == bndry) and (d != 0 or (d == 0 and parity_bit == 1)));};

      const int parity_bit = parity == FieldParity::EvenFieldParity ? (site_coords[1] & 1) : 1 - (site_coords[1] & 1);
      //
      const int my_parity    = parity == FieldParity::EvenFieldParity ? 0 : 1;
      const int other_parity = 1 - my_parity;

      Spinor res; 

      constexpr std::array<typename ArgTp::gauge_data_tp, nDir> bndr_factor{ArgTp::gauge_data_tp(1.0, 0.0),ArgTp::gauge_data_tp(-1.0, 0.0)}; 
      //
      std::array X{site_coords};	 	      
#pragma unroll
      for (int d = 0; d < nDir; d++) {
      
        const int Xd = X[d];
	// Fwd gather:
	{  

	  const bool do_halo = is_local_boundary(d, X[d], (in_accessor.extent(d) - 1), parity_bit); 
          
	  if ( do_halo ) {
	    //	
            const Link U_    = bndr_factor[d]*U(X,d, my_parity);

	    X[d] = 0;

            const Spinor in_ = in(X);
	    //
            res += U_*in_.project<+1>(d);		                  
	  } else {
            const Link U_    = U(X,d, my_parity);

	    X[d] = X[d] + (d == 0 ? parity_bit : 1);

            const Spinor in_ = in(X);
	    //		  
            res += U_*in_.project<+1>(d);		  
	  }	  
          //
          X[d] = Xd;	  
	}
	// Bwd neighbour contribution:
	{	
	  const bool do_galo = is_local_boundary(d, X[d], 0, (1-parity_bit));

          if ( do_galo ) {
            //  
	    X[d] = (in_accessor.extent(d)-1);	  

	    const Link U_    = bndr_factor[d]*U(X, d, other_parity);
	    const Spinor in_ = in(X);
            //
            res += conj(U_)*in_.project<-1>(d);              	    
          } else {  		
	    
	    X[d] = X[d] - (d == 0 ? (1- parity_bit) : 1);

	    const Link U_    = U(X,d, other_parity);
	    const Spinor in_ = in(X);
            //
	    res += conj(U_)*in_.project<-1>(d);	 
	  }
          //
          X[d] = Xd;	            
	}
      }

      return res;
    }     
    
    template<GenericSpinorFieldViewTp generic_spinor_field_view>
    void apply(auto &&transformer,
               generic_spinor_field_view &out_spinor,
               generic_spinor_field_view &in_spinor,
               generic_spinor_field_view &accum_spinor,               
               const auto cartesian_coords,
               const FieldParity parity) {	    
      // Take into account only internal points:
      // Dslash_nm = (M + 2r) \delta_nm - 0.5 * \sum_\mu  ((r - \gamma_\mu)*U_(x){\mu}*\delta_{m,n+\mu} + (r + \gamma_\mu)U^*(x-mu)_{\mu}\delta_{m,n-\mu})
      //
      // gamma_{1/2} -> sigma_{1/2}, gamma_{5} -> sigma_{3}
      //

      auto [y, x] = cartesian_coords;

      // Define accessors:
      constexpr bool is_constant = true;
      const auto U  = args.gauge.template Accessor<is_constant>();      
             
#pragma unroll
      for ( int i = 0; i < out_spinor.size(); i++ ){  	      
        auto out         = out_spinor[i].ParityAccessor();
        const auto in    = in_spinor[i].template ParityAccessor<is_constant>();
        const auto accum = accum_spinor[i].template ParityAccessor<is_constant>();        

        auto tmp = compute_parity_site_stencil(in, U, parity, {x,y});      
#pragma unroll
        for (int s = 0; s < nSpin; s++){
          out(x,y,s) = transformer(accum(x,y,s), tmp(s));//FIXME : works only for bSize = 1
        }
      }//end of for loop
    }    
    
};




