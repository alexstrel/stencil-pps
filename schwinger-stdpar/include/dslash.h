#pragma once
#include <algorithm>
#include <execution>
//
#include <dslash_helpers.h>
#include <cartesian_product.hpp>

template<typename T>
class DslashParam{
  public:
    const T M;
    const T r;    
};

template <GaugeFieldViewTp gauge_tp, int nSpin_ = 2>
class DslashArgs{
  public:
    using gauge_data_tp  = typename gauge_tp::data_tp;	  

    static constexpr std::size_t nDir   = gauge_tp::nDir;
    static constexpr std::size_t nColor = gauge_tp::nColor;
    static constexpr std::size_t nSpin  = nSpin_; 

    const gauge_tp  gauge;
    
    DslashArgs( const gauge_tp &gauge) : gauge(gauge) {}
};

template <typename Arg>
class Dslash{
  public:

    const Arg &args;

    Dslash(const Arg &args) : args(args) {}     

    template<int sign>
    inline decltype(auto) proj(const auto &in, const int dir){

      using Spinor = typename std::remove_cvref_t<decltype(in)>;
      using DataTp = typename std::remove_cvref_t<Arg>::gauge_data_tp;
      
      Spinor res;	    

      auto ic = [](auto c){ return DataTp(-c.imag(), c.real());};

      if constexpr (sign == +1) {
       switch (dir) {
          case 0 :
            res[0] = in[0] - in[1];//x+1
            res[1] = in[1] - in[0];//x+1

            break;

          case 1 :
            res[0] = in[0] + ic(in[1]);//y+1
            res[1] = in[1] - ic(in[0]);//y+1

            break;
        }
      } else if constexpr (sign == -1) {	      
        switch (dir) {
          case 0 :
            res[0] = in[0] + in[1];//x-1
            res[1] = in[1] + in[0];//x-1

            break;

          case 1 :
            res[0] = in[0] - ic(in[1]);//y-1
            res[1] = in[1] + ic(in[0]);//y-1

            break;
        }	      
      }

      return res;
    }  
    

    template<int nDir, int nSpin>
    inline decltype(auto) compute_site_stencil(auto &out, const auto &in, const auto &U, const std::array<int, nDir> site_coords){
      using ArgTp = typename std::remove_cvref_t<Arg>;

      using DataTp = ArgTp::gauge_data_tp;
      //
      using Link   = DataTp; 
      using Spinor = std::array<DataTp, nSpin>;    
    
      Spinor res;
      
      constexpr std::array<DataTp, nDir> bndr_factor{DataTp(1.0),DataTp(-1.0)}; 
      
      for (int d = 0; d < nDir; d++) {

	// Fwd gather:
	{ 
          std::array X{site_coords};	  	
          //
	  const bool ghost = false;
	  //
	  if ( ghost ) {
	    //	  
	  } else {
            const bool local_bndr =  X[d] == (in.extent(d)-1);

	    X[d] = local_bndr ? 0 : X[d] + 1;

            const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};
	    //
            const Link U_ = local_bndr ? bndr_factor[d]*U(X[0],X[1],d) : U(X[0],X[1],d);

            res += U_*proj<+1>(in_, d);		  
	  }	  
	}
	// Bwd neighbour contribution:
	{
          std::array X{site_coords};	  	
	  //
	  const bool ghost = false;
          //	
          if ( ghost ) {
            //    
          } else {  	
	    //
	    const bool local_bndr = X[d] == 0;

            X[d] = local_bndr ? (in.extent(d)-1) : X[d] - 1;		  

	    const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};

	    const Link U_ = local_bndr ? bndr_factor[d]*U(X[0],X[1],d) : U(X[0],X[1],d);

	    res += conj(U_)*proj<-1>(in_, d);	 
	  }
	}
      }               

      return res;
    }     

    template<BlockSpinorFieldViewTp  block_spinor_field_view>
    void apply(auto &&transformer,
               block_spinor_field_view &out_block_spinor,
               block_spinor_field_view &in_block_spinor,
               const auto cartesian_coords) {	    
      // Take into account only internal points:
      // Dslash_nm = (M + 2r) \delta_nm - 0.5 * \sum_\mu  ((r - \gamma_\mu)*U_(x){\mu}*\delta_{m,n+\mu} + (r + \gamma_\mu)U^*(x-mu)_{\mu}\delta_{m,n-\mu})
      //
      // gamma_{1/2} -> sigma_{1/2}, gamma_{5} -> sigma_{3}
      //
      using ArgTp = typename std::remove_cvref_t<Arg>;

      constexpr auto nSpin = ArgTp::nSpin;      
      constexpr auto nDir  = ArgTp::nDir;       

      auto [y, x] = cartesian_coords;

      // Define accessors:
      constexpr bool is_constant = true;      
      const auto U  = args.gauge.template Accessor<is_constant>();       
      
#pragma unroll
      for ( int i = 0; i < out_block_spinor.size(); i++ ){  	      
        auto out      = out_block_spinor[i].Accessor();
        const auto in = in_block_spinor[i].template Accessor<is_constant>();

        auto tmp = compute_site_stencil<nDir, nSpin>(out, in, U, {x,y});      
#pragma unroll
        for (int s = 0; s < nSpin; s++){
          out(x,y,s) = transformer(in(x,y,s), tmp[s]);
        }
      }//end of for loop
    }
    
    template<SpinorFieldViewTp spinor_field_view>    
    void apply(auto &&transformer,
	       spinor_field_view &out_spinor, 
	       spinor_field_view &in_spinor, 
	       const auto cartesian_coords, 
	       const FieldParity parity) {

      // Take into account only internal points:
      // Dslash_nm = (M + 4r) \delta_nm - 0.5 * \sum_\mu  ((r - \gamma_\mu)*U_(x){\mu}*\delta_{m,n+\mu} + (r + \gamma_\mu)U^*(x-mu)_{\mu}\delta_{m,n-\mu})
      //
      // gamma_{1/2} -> sigma_{1/2}, gamma_{5} -> sigma_{3}
      //
      using ArgTp  = typename std::remove_cvref_t<Arg>;
      	    
      using DataTp = typename ArgTp::gauge_data_tp;
      //
      constexpr auto nSpin = ArgTp::nSpin;
      //
      using Link   = DataTp; 
      using Spinor = std::array<DataTp, nSpin>;
      //
      auto is_local_boundary = [](const auto d, const auto coord, const auto bndry, const auto parity_bit){ 
	      return ((coord == bndry) && (d != 0 || (d == 0 && parity_bit)));};

      auto [y, x] = cartesian_coords;

      // Define accessors:
      constexpr bool is_constant = true;
      
      auto out      = out_spinor.Accessor();
      const auto in = in_spinor.template Accessor<is_constant>();
      //
      const auto U  = args.gauge.template ExtAccessor<is_constant>();
      
      const int parity_bit = parity == FieldParity::EvenFieldParity ? (y % 1) : 1 - (y % 1);
      //
      const int my_parity    = parity == FieldParity::EvenFieldParity ? 0 : 1;
      const int other_parity = 1 - my_parity;

      Spinor tmp;

      constexpr auto nDir = ArgTp::nDir; 

      constexpr std::array<DataTp, nDir> bndr_factor{DataTp(1.0),DataTp(-1.0)}; 
#pragma unroll
      for (int d = 0; d < nDir; d++) {

	// Fwd gather:
	{  
          std::array<int, nDir> X{x, y};	 	

	  const bool ghost = false;
          
	  if ( ghost ) {
	    //	  
	  } else {
	    const bool local_boundary_flag = is_local_boundary(d, X[d], (in.extent(d)-1), parity_bit); 

	    X[d] = (X[d] + (d == 0 ? parity_bit : 1) + in.extent(d)) % in.extent(d);

            const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};
	    //
            const Link U_ = local_boundary_flag ? (bndr_factor[d]*U(X[0],X[1],d, my_parity)) : U(X[0],X[1],d, my_parity);

            tmp += U_*proj<+1>(in_, d);		  
	  }	  
	}
	// Bwd neighbour contribution:
	{
          std::array<int, nDir> X{x, y};	
          //
	  const bool ghost = false;

          if ( ghost ) {
            //    
          } else {  		
            const bool local_boundary_flag = is_local_boundary(d, X[d], 0, (1-parity_bit));
	    
	    X[d] = (X[d] - (d == 0 ? (1- parity_bit) : 1) + in.extent(d)) % in.extent(d);	  

	    const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};
            //
	    const Link U_ = local_boundary_flag ? bndr_factor[d]*U(X[0],X[1],d, other_parity) : U(X[0],X[1],d, other_parity);

	    tmp += conj(U_)*proj<-1>(in_, d);	 
	  }
	}
      }

#pragma unroll
      for (int s = 0; s < nSpin; s++)
        out(x,y,s) = transformer(in(x,y,s), tmp[s]);
    }
    
};




