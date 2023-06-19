#pragma once

#include <cstdio>
#include <iostream>
#include <complex>
//
#include <matrix.h>


namespace impl{


  template <ComplexTp T, int Nc, int Ns, int bSize = 1>
    class ColorSpinor {
      public:
        static constexpr int N = Nc*Ns;
              
        using View3D   =  stdex::mdspan<T,       stdex::extents<indx_type, Nc, Ns, bSize>, stdex::layout_right, stdex::default_accessor<T>>;
        using CView3D  =  stdex::mdspan<const T, stdex::extents<indx_type, Nc, Ns, bSize>, stdex::layout_right, stdex::default_accessor<const T>>;

        using View2D   =  stdex::mdspan<T,       stdex::extents<indx_type, N, bSize>, stdex::layout_right, stdex::default_accessor<T>>;
        using CView2D  =  stdex::mdspan<const T, stdex::extents<indx_type, N, bSize>, stdex::layout_right, stdex::default_accessor<const T>>;      
        //
        std::array<T, N, bSize> data = {};

        inline ColorSpinor<T, Nc, Ns, bSize>() {
#pragma unroll
          for (int i = 0; i < data.size(); i++) { data[i] = T{0, 0}; }
        }

        inline ColorSpinor<T, Nc, Ns, bSize> operator-() const
        {
          ColorSpinor<T, Nc, Ns, bSize> a;
#pragma unroll
          for (int i = 0; i < data.size(); i++) { a.data[i] = -data[i]; }
          return a;
        }

        inline ColorSpinor<T, Nc, Ns, bSize>& operator+=(const ColorSpinor<T, Nc, Ns, bSize> &a) {
#pragma unroll
          for (int i = 0; i < data.size(); i++) { data[i] = data[i] + a.data[i]; }
          return *this;
        }

        inline ColorSpinor<T, Nc, Ns, bSize> &operator*=(const ComplexTp auto &a)
        {
#pragma unroll
          for (int i = 0; i < data.size(); i++) { data[i] = data[i]*a; }
          return *this;
        }

        inline ColorSpinor<T, Nc, Ns, bSize> &operator*=(const auto FloatTp &a)
        {
#pragma unroll
          for (int i = 0; i < data.size(); i++) { data[i] = T{data[i].real()*a, data[i].real()*a}; }
          return *this;
        }


        inline ColorSpinor<T, Nc, Ns, bSize> &operator-=(const ColorSpinor<T, Nc, Ns, bSize> &a)
        {
          if (this != &a) {
#pragma unroll
            for (int i = 0; i < data.size(); i++) { data[i] = data[i] - a.data[i]; }
          }
          return *this;
        }
      
        inline decltype(auto) cview() const {
          return CView3D(data.data());
        }

        inline decltype(auto) view() {
          return View3D(data.data());
        }

        inline decltype(auto) flat_cview() const {
          return CView2D(data.data());
        }

        inline decltype(auto) flat_view() {
          return View2D(data.data());
        }
    };


    /**
       This is the specialization for Nspin=2.  For fields with two
       spins we can define a spin projection operation.
    */
    template <ComplexTp T, int Nc, int bSize = 1> 
      class ColorSpinor<T, Nc, 2, bSize> {
      
        public: 
          static constexpr int Ns = 2;
          static constexpr int N  = Nc * Ns;
        
          using View3D   =  stdex::mdspan<T,       stdex::extents<indx_type, Nc, Ns, bSize>, stdex::layout_right, stdex::default_accessor<T>>;
          using CView3D  =  stdex::mdspan<const T, stdex::extents<indx_type, Nc, Ns, bSize>, stdex::layout_right, stdex::default_accessor<const T>>;

          using View2D   =  stdex::mdspan<T,       stdex::extents<indx_type, N, bSize>, stdex::layout_right, stdex::default_accessor<T>>;
          using CView2D  =  stdex::mdspan<const T, stdex::extents<indx_type, N, bSize>, stdex::layout_right, stdex::default_accessor<const T>>;      
      
          std::array<T, N, bSize> data = {};

          inline ColorSpinor<T, Nc, 2, bSize>()
          {
#pragma unroll
            for (int i = 0; i < data.size(); i++) { data[i] = T{0, 0}; }
          }
          
          inline ColorSpinor<T, Nc, 2, bSize>& operator+=(const ColorSpinor<T, Nc, 2, bSize> &a) {
#pragma unroll
            for (int i = 0; i < size; i++) { data[i] = data[i] + a.data[i]; }
            return *this;
          }

          inline ColorSpinor<T, Nc, 2, bSize> &operator*=(const ComplexTp auto &a)
          {
#pragma unroll
            for (int i = 0; i < size; i++) { data[i] = a*data[i]; }
            return *this;
          }          

          inline ColorSpinor<T, Nc, 2, bSize> &operator*=(const FloatTp auto &a)
          {
#pragma unroll
            for (int i = 0; i < size; i++) { data[i] = a*data[i]; }
            return *this;
          }          

          template<int sign>
          inline decltype(auto) project(const int dir){

          using data_tp = typename T::value_type; 
      
          std::array<T, N, bSize> proj;	    

          auto I = [](auto a){ return T{-a.imag(), a.real()};};

          const auto in = CView3D{data.data()}; 

          if constexpr (sign == +1) {
            switch (dir) {
              case 0 :
#pragma unroll              
              for(int c = 0; c < Nc; c++ ) {
#pragma unroll              
                for(int b = 0; b < bSize; b++) {
                  proj(c, 0, b) = in(c, 0, b) - in(c, 1, b);
                  proj(c, 1, b) = in(c, 1, b) - in(c, 0, b);
                }
              }

              break;

              case 1 :
#pragma unroll              
              for(int c = 0; c < Nc; c++ ) {
#pragma unroll              
                for(int b = 0; b < bSize; b++) {
                  proj(c, 0, b) = in(c, 0, b) + I(in(c, 1, b));
                  proj(c, 1, b) = in(c, 1, b) - I(in(c, 0, b));
                }
              }              
              break;
             }
           } else if constexpr (sign == -1) {	      
            switch (dir) {
              case 0 :
#pragma unroll              
              for(int c = 0; c < Nc; c++ ) {
#pragma unroll              
                for(int b = 0; b < bSize; b++) {
                  proj(c, 0, b) = in(c, 0, b) + in(c, 1, b);
                  proj(c, 1, b) = in(c, 1, b) + in(c, 0, b);
                }
              }

              break;

              case 1 :
#pragma unroll              
              for(int c = 0; c < Nc; c++ ) {
#pragma unroll              
                for(int b = 0; b < bSize; b++) {
                  proj(c, 0, b) = in(c, 0, b) - I(in(c, 1, b));
                  proj(c, 1, b) = in(c, 1, b) + I(in(c, 0, b));
                }
              }              
              break;
             }	      
           }
           //
           return res;
         }
  };
  
  /**
     @brief Compute the matrix-vector product y = A * x
     @param[in] A Input matrix
     @param[in] x Input vector
     @return The vector A * x
  */
  template<ComplexTp T, int Nc, int Ns, int bSize = 1>
    inline ColorSpinor<T,Nc,Ns, bSize> operator*(const Matrix<T,Nc,bSize> &A, const ColorSpinor<T,Nc,Ns,bSize> &x) {
      //
      using data_type = T::value_type;
      ColorSpinor<T,Nc,Ns,bSize> y;
      //
      auto y_ = y.view();
      //
      const auto A_ = A.cview();
      const auto x_ = x.cview();    
      //
#pragma unroll
      for (int c=0; c<Nc; c++) {
#pragma unroll
        for (int s=0; s<Ns; s++) {
          data_tp re;
          data_tp im;
#pragma unroll        
          for (int b=0; b<bSize; b++) {
	    re  = A_(c,0,b).real() * x_(0, s, b).real();
	    re -= A_(c,0,b).imag() * x_(0, s, b).imag();
	    im  = A_(c,0,b).real() * x_(0, s, b).imag();
	    im += A_(c,0,b).imag() * x_(0, s, b).real();
	    y_(c,0,b) = T{re, im};
	  }
        }
#pragma unroll
        for (int cc=1; cc<Nc; cc++) {
#pragma unroll
	  for (int s=0; s<Ns; s++) {
#pragma unroll        
            for (int b=0; b<bSize; b++) {
	      re += A_(c,cc,b).real() * x_(cc, s, b).real();
	      re -= A_(c,cc,b).imag() * x_(cc, s, b).imag();
	      im += A_(c,cc,b).real() * x_(cc, s, b).imag();
	      im += A_(c,cc,b).imag() * x_(cc, s, b).real();
	      y_(c,cc,b) += T{re, im};
	    }	  
	  }
        }
      }

      return y;
    }
  
   /**
     @brief caxpy operation on ColorSpinor objects
     @param[in] a complex scalar
     @param[in] x Vector that is scaled
     @param[in,out] y Accumulation vector
  */
  template <ComplexTp T, int Nc, int Ns, int bSize = 1>
  inline void caxpy(const T &a, const ColorSpinor<T, Nc, Ns, bSize> &x, ColorSpinor<T, Nc, Ns, bSize> &y)
  {
    auto       y_ = y.flat_view();
    const auto x_ = x.flat_cview();    
#pragma unroll
    for (int i = 0; i < Nc*Ns; i++) {
#pragma unroll
      for (int b = 0; b < bSize; b++) {    
        y_(i,b).real( a.real() * x_(i,b).real() + y_(i,b).real());
        y_(i,b).real(-a.imag() * x_(i,b).imag() + y_(i,b).real());
        y_(i,b).imag( a.imag() * x_(i,b).real() + y_(i,b).imag());
        y_(i,b).imag( a.real() * x_(i,b).imag() + y_(i,b).imag());
      }
    }
  }
  
  /**
     @brief ColorSpinor addition operator
     @param[in] x Input vector
     @param[in] y Input vector
     @return The vector x + y
  */
  template<ComplexTp T, int Nc, int Ns, int bSize = 1>
  inline ColorSpinor<T,Nc,Ns,bSize> operator+(const ColorSpinor<T,Nc,Ns,bSize> &x, const ColorSpinor<T,Nc,Ns,bSize> &y) {

    ColorSpinor<T,Nc,Ns,bSize> z;

    const auto y_ = y.flat_cview();
    const auto x_ = x.flat_cview();    

    auto z_ = z.flat_view();    

#pragma unroll
    for (int i=0; i<Nc*Ns; i++) {
#pragma unroll      
      for (int b=0; b<bSize; b++) {
	z_(i,b) = x_(i,b) + y_(i,b);
      }
    }

    return z;
  }

  /**
     @brief ColorSpinor subtraction operator
     @param[in] x Input vector
     @param[in] y Input vector
     @return The vector x + y
  */
  template<ComplexTp T, int Nc, int Ns, int bSize = 1>
    inline ColorSpinor<T,Nc,Ns,bSize> operator-(const ColorSpinor<T,Nc,Ns,bSize> &x, const ColorSpinor<T,Nc,Ns,bSize> &y) {

    ColorSpinor<T,Nc,Ns,bSize> z;

    const auto y_ = y.flat_cview();
    const auto x_ = x.flat_cview();    

    auto z_ = z.flat_view();    

#pragma unroll
    for (int i=0; i<Nc*Ns; i++) {
#pragma unroll      
      for (int b=0; b<bSize; b++) {
	z_(i,b) = x_(i,b) - y_(i,b);
      }
    }

    return z;
  }

  /**
     @brief Compute the scalar-vector product y = a * x
     @param[in] a Input scalar
     @param[in] x Input vector
     @return The vector a * x
  */
  template<ArithmeticTp S, ComplexTp T, int Nc, int Ns, int bSize = 1> 
    inline ColorSpinor<T,Nc,Ns,bSize> operator*(const S &a, const ColorSpinor<T,Nc,Ns,bSize> &x) {

    ColorSpinor<T,Nc,Ns,bSize> y;

    auto y_ = y.flat_cview();
    const auto x_ = x.flat_cview();    

#pragma unroll
    for (int i=0; i<Nc*Ns; i++) {
#pragma unroll      
      for (int b=0; b<bSize; b++) {
	y_(i,b) = a*x_(i,b);
      }
    }

    return y;
  }
  
   
} // 
