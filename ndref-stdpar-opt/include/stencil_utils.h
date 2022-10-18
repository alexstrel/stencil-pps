#pragma once
#include <common.h>
#include <field.h>

constexpr long double pi{std::numbers::pi_v<long double>};
//constexpr long double _2pi = (2*pi);

template<typename Policy,  ArithmeticTp T, bool is_zero, int... M>
void create_field(const Policy &p,
                  std::vector<impl::StencilCell<T, M...>> &out,
                  const std::array<int, 3> nl,
                  const T kappa,
                  const T length,
                  const T time) {
                  
  std::vector<impl::StencilCell<T, M...>> tmp(out.size());                  

  printf("Stencil size %d\n", tmp[0].size());
  
  const std::array<T, 3> dnl{length / (nl[0]+static_cast<T>(1.0)), length /(nl[1]+static_cast<T>(1.0)), length / (nl[2]+static_cast<T>(1.0))} ; 	
  
  const T exponent  = is_zero ? 0.0 : exp(-3.0*kappa*(pi*pi)*time/(length*length));
  
  T z = dnl[2];
  if constexpr (sizeof...(M) != nl.size()) {
    std::cout << "Dimension is not supported" << std::endl;
  }
  constexpr std::array<int, sizeof...(M)> cell_vol{M...};
  std::array<int, sizeof...(M)> domn_vol;  
  
  for(int i = 0; i < nl.size(); i++){
    domn_vol[i] = nl[i] / cell_vol[i];
  }
  
  for (int k = 0; k < nl[2]; k++){
    T y = dnl[1];
    int kg = k % domn_vol[2];
    int kl = k / domn_vol[2];
    for (int j = 0; j < nl[1]; j++){
      T x = dnl[0];
      int jg = j % domn_vol[1];
      int jl = j / domn_vol[1];
      for (int i = 0; i < nl[0]; i++){
        int ig = i % domn_vol[0];
        int il = i / domn_vol[0];        
        //
        int s = k*nl[0]*nl[1] + j*nl[0] + i;
        T val = exponent * sin(pi*x/length) * sin(pi*y/length)* sin(pi*z/length);
        //
        int sg = kg*domn_vol[0]*domn_vol[1] + jg*domn_vol[0] + ig;
        int sl = kl*cell_vol[0]*cell_vol[1] + jl*cell_vol[0] + il;
        tmp[sg][sl] = val;
        x += dnl[0];
      }
      y += dnl[1];
    }
    z += dnl[2];
  }
  
  std::copy(p, tmp.begin(), tmp.end(), out.begin());
  
  return;
}

template<typename Policy,  ArithmeticTp T>
void create_field_v2(const Policy &p,
                     std::vector<T> &out,
                     const std::array<int, 3> nl,
                     const T kappa,
                     const T length,
                     const T time) {
                  
  std::vector<T> tmp(out.size());                  
  
  const std::array<T, 3> dnl{length / (nl[0]+static_cast<T>(1.0)), length /(nl[1]+static_cast<T>(1.0)), length / (nl[2]+static_cast<T>(1.0))} ; 	
  
  const T exponent  = exp(-3.0*kappa*(pi*pi)*time/(length*length));
  
  T z = dnl[2];
  
  for (int k = 0; k < nl[2]; k++){
    T y = dnl[1];
    for (int j = 0; j < nl[1]; j++){
      T x = dnl[0];
      for (int i = 0; i < nl[0]; i++){
        int s = k*nl[0]*nl[1] + j*nl[0] + i;
        tmp[s] = exponent * sin(pi*x/length) * sin(pi*y/length)* sin(pi*z/length);
        x += dnl[0];
      }
      y += dnl[1];
    }
    z += dnl[2];
  }
  
  std::copy(p, tmp.begin(), tmp.end(), out.begin());
  
  return;
}

template<ArithmeticTp T, bool rel_error, int... M>
T accuracy(const std::vector<impl::StencilCell<T, M...>> &ref, std::vector<impl::StencilCell<T, M...>> &sol) {

  double err       = 0.0;
  double ref_norm  = 0.0;
  double sol_norm  = 0.0;

  for (int i = 0; i < ref.size(); i++){
    for(int j = 0; j < ref[0].size(); j++) {
      err      += (ref[i][j] - sol[i][j])*(ref[i][j] - sol[i][j]);
      ref_norm += ref[i][j]*ref[i][j];
      sol_norm += sol[i][j]*sol[i][j];
    }
  }
 
  std::cout <<  std::setprecision(15) << "Solution norm" << sol_norm <<", Ref norm " << ref_norm << std::endl;
  
  if constexpr (rel_error) return (T)sqrt(err/ref_norm);
  else                     return (T)sqrt(err);
}

