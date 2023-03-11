#pragma once
#include <algorithm>
#include <execution>
//
#include <dslash_helpers.h>
#include <cartesian_product.hpp>

template<typename T>
class DslashParam{
  public:
    const T kappa;	
};

template <typename gauge_tp, typename param_tp, int nSpin_ = 2>
class DslashArgs{
  public:
    using gauge_data_tp  = typename gauge_tp::data_tp;	  

    static constexpr std::size_t nDir   = gauge_tp::nDir;
    static constexpr std::size_t nColor = gauge_tp::nColor;
    static constexpr std::size_t nSpin  = nSpin_; 

    const gauge_tp  gauge;
    
    const param_tp param;

    DslashArgs( const gauge_tp &gauge, const param_tp &param) : gauge(gauge), param(param) {}
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

    void apply(auto &out_spinor, const auto &in_spinor, const auto cartesian_coords) {

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
      const auto U  = args.gauge.Accessor();

      const auto kappa = args.param.kappa;

      std::array<DataTp, nSpin> tmp;

      constexpr auto nDir = std::remove_cvref_t<Arg>::nDir; 

      constexpr std::array<DataTp, nDir> bndr_factor{DataTp(1.0),DataTp(-1.0)}; 
#pragma unroll
      for (int d = 0; d < nDir; d++) {

        std::array<int, nDir> X{x, y};

	// Fwd gather:
	{   	
	  if ( X[d] == (in.extent(d)-1) ) {
	    //	  
	    X[d] = 0;

	    const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};

	    tmp += (bndr_factor[d]*U(x,y,d))*proj<+1>(in_, d);
	  } else {
	    X[d] += 1;

            const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};

            tmp += U(x,y,d)*proj<+1>(in_, d);		  
	  }	  
	}
	// Bwd neighbour contribution:
	{
          if ( X[d] == 0 ) {
            //    
	    X[d] = (in.extent(d)-1);

            const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};

            tmp += conj(bndr_factor[d]*U(X[0],X[1],d))*proj<-1>(in_, d);
          } else {  		
            X[d] -= 1;		  

	    const Spinor in_{in(X[0],X[1],0), in(X[0],X[1],1)};

	    tmp += conj(U(X[0],X[1],d))*proj<-1>(in_, d);	 
	  }
	}
      }

#pragma unroll
      for (int s = 0; s < nSpin; s++)
        out(x,y,s) = in(x,y,s) - kappa*tmp[s];
    }
};

#if 0
template<typename Kernel, typename KernelArgs>
class Mat{
  private:
    std::unique_ptr<Kernel> dslash_kernel_ptr;	

  public:

    Mat(const KernelArgs &args) : dslash_kernel_ptr(new Kernel(args)) {}

    void operator()(auto &out, auto &in){
      assert(in.GetFieldOrder() == FieldOrder::LexFieldOrder);	    
      // Take into account only internal points:
      const auto [Nx, Ny] = in.GetCBDims(); //Get CB dimensions

      auto X = std::views::iota(0, Nx-1);
      auto Y = std::views::iota(0, Ny-1);

      auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest	    
       	    
      auto DslashKernel = [&dslash_kernel = *dslash_kernel_ptr.get(), out_ = out.Get(), in_ = in.Get()] (const auto i) { dslash_kernel.apply(out_, in_, i); };    
      //
      std::for_each(std::execution::par_unseq,
                    idx.begin(),
                    idx.end(),
                    DslashKernel);       
    }
};

#endif
