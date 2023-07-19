#pragma once
#include <algorithm>
#include <execution>
#include <ranges>
//
#include <core/cartesian_product.hpp>
//
#include <fields/field_accessor.h>

template<typename T>
class DslashParam{
  public:
    const T M;
    const T r;
};

constexpr bool is_constant = true;

//!template <GaugeFieldViewTp gauge_tp, int nSpin_ = 2, int bSize_  = 1>
template <GaugeFieldViewTp gauge_tp, int nSpin_ = 2>
class DslashArgs {
  public:
    using gauge_data_tp  = typename gauge_tp::data_tp;	  

    static constexpr std::size_t nDir   = gauge_tp::Ndir();
    static constexpr std::size_t nDim   = gauge_tp::Ndim();    
    
//!    static constexpr std::size_t bSize  = bSize_;    

    static constexpr std::size_t bSize  = 1;    
    
    using LinkAccessor = FieldAccessor<gauge_tp, is_constant, bSize>;//only constant links   
    using LinkTp       = LinkAccessor::LinkTp;

    const LinkAccessor U;

    DslashArgs( const gauge_tp &gauge) : U(gauge) {}
};


template <typename Arg>
class Dslash{
  public:
    using ArgTp  = typename std::remove_cvref_t<Arg>;

    const Arg &args;

    Dslash(const Arg &args) : args(args) {}     

    /** Convert tuple into std::array in the inverse order:
    */
    template<std::size_t... Id>
    inline decltype(auto) get_cartesian_coords(std::index_sequence<Id...>, const auto& x) const {
      return std::array<int, sizeof... (Id)>{{std::get<sizeof...(Id)-Id-1>(x)... }};
    }

    inline decltype(auto) convert_coords(const auto &x) const {
      return get_cartesian_coords(std::make_index_sequence<ArgTp::nDim>{}, x);
    }

    template<bool dagger>
    inline decltype(auto) compute_parity_site_stencil(const auto &in, const FieldParity parity, std::array<int, ArgTp::nDim> site_coords){
    
      using Link   = ArgTp::LinkTp; 
      using Spinor = typename std::remove_cvref_t<decltype(in)>::SpinorTp;
      
      //Define accessor wrappers:    
      auto is_local_boundary = [](const auto d, const auto coord, const auto bndry, const auto parity_bit){ 
	      return ((coord == bndry) and (d != 0 or (d == 0 and parity_bit == 1)));};

      const int parity_bit = parity == FieldParity::EvenFieldParity ? (site_coords[1] & 1) : 1 - (site_coords[1] & 1);
      //
      const int my_parity    = parity == FieldParity::EvenFieldParity ? 0 : 1;
      const int other_parity = 1 - my_parity;

      Spinor res; 

      constexpr std::array<typename ArgTp::gauge_data_tp, ArgTp::nDir> bndr_factor{ArgTp::gauge_data_tp(1.0, 0.0),ArgTp::gauge_data_tp(-1.0, 0.0)}; 
      //
      auto X = site_coords | std::views::all;

#pragma unroll
      for (int d = 0; d < ArgTp::nDir; d++) {
      
        const int Xd = X[d];
	// Fwd gather:
	{  

          constexpr int sign = dagger ? -1 : +1;

	  const bool do_halo = is_local_boundary(d, X[d], (in.Extent(d) - 1), parity_bit); 
          
	  if ( do_halo ) {
	    //	
            const Link U_    = bndr_factor[d]*args.U(X,d, my_parity);

	    X[d] = 0;

            const Spinor in_ = in(X);
	    //
            res += U_*in_.project<sign>(d);		                  
	  } else {
            const Link U_    = args.U(X,d, my_parity);

	    X[d] = X[d] + (d == 0 ? parity_bit : 1);

            const Spinor in_ = in(X);
	    //		  
            res += U_*in_.project<sign>(d);		  
	  }	  
          //
          X[d] = Xd;	  
	}
	// Bwd neighbour contribution:
	{
	  const bool do_halo = is_local_boundary(d, X[d], 0, (1-parity_bit));

          constexpr int sign = dagger ? +1 : -1;

          if ( do_halo ) {
            //  
	    X[d] = (in.Extent(d)-1);	  

	    const Link U_    = bndr_factor[d]*args.U(X, d, other_parity);
	    const Spinor in_ = in(X);
            //
            res += conj(U_)*in_.project<sign>(d);              	    
          } else {  		
	    
	    X[d] = X[d] - (d == 0 ? (1- parity_bit) : 1);

	    const Link U_    = args.U(X,d, other_parity);
	    const Spinor in_ = in(X);
            //
	    res += conj(U_)*in_.project<sign>(d);	 
	  }
          //
          X[d] = Xd;	            
	}
      }

      return res;
    }     
 
    template <bool dagger>   
    void apply(GenericSpinorFieldViewTp auto &out_spinor,
               const GenericSpinorFieldViewTp auto &in_spinor,
               const GenericSpinorFieldViewTp auto &aux_spinor,
               auto &&post_transformer,               
               const auto cartesian_idx,
               const FieldParity parity) {	    
      // Dslash_nm = (M + 2r) \delta_nm - 0.5 * \sum_\mu  ((r - \gamma_\mu)*U_(x){\mu}*\delta_{m,n+\mu} + (r + \gamma_\mu)U^*(x-mu)_{\mu}\delta_{m,n-\mu})
      //
      // gamma_{1/2} -> sigma_{1/2}, gamma_{5} -> sigma_{3}
      //
      using S = typename std::remove_cvref_t<decltype(out_spinor[0])>;       

      auto X  = convert_coords(cartesian_idx);

      auto X_view = X | std::views::all; 

#pragma unroll
      for ( int i = 0; i < out_spinor.size(); i++ ){  	      
      
        auto out         = FieldAccessor<S>{out_spinor[i]};
        const auto in    = FieldAccessor<S, is_constant>{in_spinor[i]};
        const auto aux   = FieldAccessor<S, is_constant>{aux_spinor[i]};        
        //
        auto res = compute_parity_site_stencil<dagger>(in, parity, X);
        //
        const auto aux_  = aux(X_view);
        //
        post_transformer(aux_, res);
        //
#pragma unroll
        for (int s = 0; s < S::Nspin(); s++){
          out(X_view, s) = res(s);
        }        
      }//end of for loop
    }    

    template <bool dagger>   
    void apply(GenericSpinorFieldViewTp auto &out_spinor,
               const GenericSpinorFieldViewTp auto &in_spinor,
               const auto cartesian_idx,
               const FieldParity parity) {	    
      // Dslash_nm = \sum_\mu  ((r - \gamma_\mu)*U_(x){\mu}*\delta_{m,n+\mu} + (r + \gamma_\mu)U^*(x-mu)_{\mu}\delta_{m,n-\mu})
      //
      // gamma_{1/2} -> sigma_{1/2}, gamma_{5} -> sigma_{3}
      //
      using S = typename std::remove_cvref_t<decltype(out_spinor[0])>; 

      auto X  = convert_coords(cartesian_idx);

      auto X_view = X | std::views::all;

#pragma unroll
      for ( int i = 0; i < out_spinor.size(); i++ ){  	      
      
        auto out         = FieldAccessor<S>{out_spinor[i]};
        const auto in    = FieldAccessor<S, is_constant>{in_spinor[i]};
        //
        auto res = compute_parity_site_stencil<dagger>(in, parity, X);
    
#pragma unroll
        for (int s = 0; s < S::Nspin(); s++){
          out(X_view, s) = res(s);//FIXME : works only for bSize = 1
        }
      }//end of for loop
    }    

    
};




