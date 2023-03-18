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

template <typename gauge_tp, int nSpin_ = 2>
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

    void apply(auto &&transformer, 
	       auto &out_spinor, 
	       const auto &in_spinor, 
	       const auto cartesian_coords) {
      // Take into account only internal points:
      // Dslash_nm = (M + 2r) \delta_nm - 0.5 * \sum_\mu  ((r - \gamma_\mu)*U_(x){\mu}*\delta_{m,n+\mu} + (r + \gamma_\mu)U^*(x-mu)_{\mu}\delta_{m,n-\mu})
      //
      // gamma_{1/2} -> sigma_{1/2}, gamma_{5} -> sigma_{3}
      //

      using ArgTp = typename std::remove_cvref_t<Arg>;

      using DataTp = ArgTp::gauge_data_tp;
      //
      constexpr auto nSpin = ArgTp::nSpin;
      //
      using Link   = DataTp; 
      using Spinor = std::array<DataTp, nSpin>;

      auto [y, x] = cartesian_coords;

      // Define accessors:
      auto out      = out_spinor.Accessor();
      const auto in = in_spinor.Accessor();
      const auto U  = args.gauge.Accessor();

      Spinor tmp;

      constexpr auto nDir = ArgTp::nDir; 

      constexpr std::array<DataTp, nDir> bndr_factor{DataTp(1.0),DataTp(-1.0)}; 

#pragma unroll
      for (int d = 0; d < nDir; d++) {

	// Fwd gather:
	{ 
          std::array<int, nDir> X{x, y};	  	
          //
	  if ( X[d] == (in.extent(d)-1) ) {
	    //	  
	    X[d] = 0;

	    const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};

	    const Link U_ = U(x,y,d);

	    tmp += (bndr_factor[d]*U_)*proj<+1>(in_, d);

	  } else {
	    X[d] += 1;

            const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};

	    const Link U_ = U(x,y,d);

            tmp += U_*proj<+1>(in_, d);		  
	  }	  
	}
	// Bwd neighbour contribution:
	{
          std::array<int, nDir> X{x, y};	  	
          //	
          if ( X[d] == 0 ) {
            //    
	    X[d] = (in.extent(d)-1);

            const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};

	    const Link U_ = U(X[0],X[1],d);

            tmp += conj(bndr_factor[d]*U_)*proj<-1>(in_, d);
          } else {  	
	    //	  
            X[d] -= 1;		  

	    const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};

	    const Link U_ = U(X[0],X[1],d);

	    tmp += conj(U_)*proj<-1>(in_, d);	 
	  }
	}
      }      
      
#pragma unroll
      for (int s = 0; s < nSpin; s++){
        out(x,y,s) = transformer(in(x,y,s), tmp[s]);
      }
    }
    
    template<bool do_pre_transform = false>
    void apply(auto &&pre_transformer,
	       auto &&post_transformer, 
	       auto &out_spinor, 
	       const auto &in_spinor, 
	       const auto cartesian_coords, 
	       const FieldParity parity) {

      // Take into account only internal points:
      // Dslash_nm = (M + 4r) \delta_nm - 0.5 * \sum_\mu  ((r - \gamma_\mu)*U_(x){\mu}*\delta_{m,n+\mu} + (r + \gamma_\mu)U^*(x-mu)_{\mu}\delta_{m,n-\mu})
      //
      // gamma_{1/2} -> sigma_{1/2}, gamma_{5} -> sigma_{3}
      //
	    
      using DataTp = typename std::remove_cvref_t<Arg>::gauge_data_tp;
      //
      constexpr auto nSpin = std::remove_cvref_t<Arg>::nSpin;
      //
      using Link   = DataTp; 
      using Spinor = std::array<DataTp, nSpin>;

      auto [y, x] = cartesian_coords;

      // Define accessors:
      auto out      = out_spinor.Accessor();
      const auto in = in_spinor.Accessor();
      const auto U  = args.gauge.ExtAccessor();

      const auto kappa = 1.0;//
      
      const int parity_bit = parity == FieldParity::EvenFieldParity ? (y % 1) : 1 - (y % 1);
      //
      const int my_parity    = parity == FieldParity::EvenFieldParity ? 0 : 1;
      const int other_parity = 1 - my_parity;

      std::array<DataTp, nSpin> tmp;

      constexpr auto nDir = std::remove_cvref_t<Arg>::nDir; 

      constexpr std::array<DataTp, nDir> bndr_factor{DataTp(1.0),DataTp(-1.0)}; 
#pragma unroll
      for (int d = 0; d < nDir; d++) {
      
        const bool fwd_ghost_flag = d != 0 || (d == 0 && parity_bit);      
        const bool bwd_ghost_flag = d != 0 || (d == 0 && (1-parity_bit));              

        std::array<int, nDir> X{x, y};               

	// Fwd gather:
	{   	
	  if ( X[d] == (in.extent(d)-1) && fwd_ghost_flag) {
	    //	  
	    X[d] = 0;

	    const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};

	    tmp += (bndr_factor[d]*U(x,y,d,my_parity))*proj<+1>(in_, d);
	  } else {
	    X[d] += (d == 0 ? parity_bit : 1);

            const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};

            tmp += U(x,y,d, my_parity)*proj<+1>(in_, d);		  
	  }	  
	}
	// Bwd neighbour contribution:
	{
          if ( X[d] == 0 && bwd_ghost_flag) {
            //    
	    X[d] = (in.extent(d)-1);

            const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};

            tmp += conj(bndr_factor[d]*U(X[0],X[1],d))*proj<-1>(in_, d);
          } else {  		
            X[d] -= (d == 0 ? (1- parity_bit) : 1);		  

	    const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};

	    tmp += conj(U(X[0],X[1],d, other_parity))*proj<-1>(in_, d);	 
	  }
	}
      }

#pragma unroll
      for (int s = 0; s < nSpin; s++)
        out(x,y,s) = in(x,y,s) - kappa*tmp[s];
    }
    
};




