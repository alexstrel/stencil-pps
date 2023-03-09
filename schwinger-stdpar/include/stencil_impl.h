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

template <typename gauge_tp, typename param_tp >
class StencilArgs{
  public:
    using gauge_data_tp  = typename gauge_tp::data_tp;	  

    const gauge_tp  gauge;
    
    const param_tp param;

    StencilArgs( const gauge_tp &gauge, const param_tp &param) : gauge(gauge), param(param) {}
};

template< ComplexTp T>
inline decltype(auto) operator*(const T &scalar, const std::array<T,2> &a){
  std::array<T,2> result;
#pragma unroll
  for(int i = 0; i < 2; i++){
    result[i] = scalar.real() * a[i].real() - scalar.imag() * a[i].imag();
    result[i] = scalar.real() * a[i].imag() + scalar.imag() * a[i].real();
  }

  return result;
}

template< ComplexTp T>
inline decltype(auto) operator+=(std::array<T,2> &a, const std::array<T,2> &b){
 
#pragma unroll
  for(int i = 0; i < 2; i++){
    a[i] = a[i].real() * b[i].real() - a[i].imag() * b[i].imag();
    a[i] = a[i].real() * b[i].imag() + a[i].imag() * b[i].real();
  }

  return a;
}


template <typename Arg>
class Stencil{
  public:

    const Arg &args;

    Stencil(const Arg &args) : args(args) {}     

    void apply(auto &out, const auto &in, const auto cartesian_coords) {

      // Take into account only internal points:
      // Dslash_nm = (M + 4r) \delta_nm - 0.5 * \sum_\mu  ((r - \gamma_\mu)*U_(x){\mu}*\delta_{m,n+\mu} + (r + \gamma_\mu)U^*(x-mu)_{\mu}\delta_{m,n-\mu})
      //
      // gamma_{1/2} -> sigma_{1/2}, gamma_{5} -> sigma_{3}
      //
	    
      using DataTp = typename std::remove_cvref_t<Arg>::gauge_data_tp;
      using Link   = DataTp; 
      using Spinor = std::array<DataTp, 2>;

      constexpr int spin_0 = 0;
      constexpr int spin_1 = 1;

      auto ix = [=](auto x){ 
        return DataTp(-x.imag(), x.real());
      };

      auto [y, x] = cartesian_coords;
      
      auto proj = [=](const auto& s, const int dir, const int sign) {
         Spinor res;

	 switch (dir) {
	   case 0 :
	     switch (sign) {
	       case +1 : 
	         res[0] = s(x+1,y,spin_0) - s(x+1,y,spin_1);
		 res[1] = s(x+1,y,spin_1) - s(x+1,y,spin_0);
                 break;
               case -1 :
                 res[0] = s(x-1,y,spin_0) + s(x-1,y,spin_1); 		 
		 res[1] = s(x-1,y,spin_1) + s(x-1,y,spin_0);
		 break;
	     }

           case 1 :
             switch (sign) {
               case +1 :
                 res[0] = s(x,y+1,spin_0) + ix(s(x,y+1,spin_1));
		 res[1] = s(x,y+1,spin_1) - ix(s(x,y+1,spin_0));
                 break;
               case -1 :
                 res[0] = s(x,y-1,spin_0) - ix(s(x,y-1,spin_1));
		 res[1] = s(x,y-1,spin_1) + ix(s(x,y-1,spin_0));
                 break;
             }              	     
	 }
	       
	 return res;
      };
 
      auto o       = out.Accessor();
      const auto i = in.Accessor();      

      const auto u = args.gauge.Accessor();

      const auto kappa = args.param.kappa;

      std::array<DataTp, 2> s;

      constexpr int nDir = 2; 
#pragma unroll
      for (int d = 0; d < nDir; d++) {
	// Forward hop:      
	s += u(x,y,d)*proj(i,d,+1);

	// Backward hop:
        s += conj(u(x-1,y,d))*proj(i,d,-1);	 
      }

      o(x,y,spin_0) = i(x,y,spin_0) - kappa*s[spin_0];
      o(x,y,spin_1) = i(x,y,spin_1) - kappa*s[spin_1];
    }
};


template<typename Kernel, typename KernelArgs>
class Mat{
  private:
    std::unique_ptr<Kernel> stencil_kernel_ptr;	

  public:

    Mat(const KernelArgs &args) : stencil_kernel_ptr(new Kernel(args)) {}

    void operator()(auto &out, auto &in){
      assert(in.GetFieldOrder() == FieldOrder::LexFieldOrder);	    
      // Take into account only internal points:
      const auto [Nx, Ny] = in.GetCBDims(); //Get CB dimensions

      auto X = std::views::iota(1, Nx-2);
      auto Y = std::views::iota(1, Ny-2);

      auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest	    
       	    
      auto StencilKernel = [&stencil_kernel = *stencil_kernel_ptr.get(), out_ = out.Get(), in_ = in.Get()] (const auto i) { stencil_kernel.apply(out_, in_, i); };    
      //
      std::for_each(std::execution::par_unseq,
                    idx.begin(),
                    idx.end(),
                    StencilKernel);       
    }
};


