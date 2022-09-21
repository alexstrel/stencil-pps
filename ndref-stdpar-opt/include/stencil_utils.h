#pragma once
#include <common.h>
#include <field.h>

const double _pi  = acos(-1);
const double _2pi = (2*_pi);

template<typename Policy,  ArithmeticTp T, int... M>
void create_field(const Policy &p,
                  std::vector<impl::StencilCell<T, M...>> &out,
                  const std::array<int, 3> nl,
                  const T kappa,
                  const T length,
                  const T time) {
                  
  std::vector<impl::StencilCell<T, M...>> tmp(out.size());                  
  
  const std::array<T, 3> dnl{length / (nl[0]+static_cast<T>(1.0)), length /(nl[1]+static_cast<T>(1.0)), length / (nl[2]+static_cast<T>(1.0))} ; 	
  
  const T exponent  = exp(-3.0*kappa*_pi*_pi*time/(length*length));
  
  T z = dnl[2];
  
  for (int k = 0; k < nl[2]; k++){
    T y = dnl[1];
    for (int j = 0; j < nl[1]; j++){
      T x = dnl[0];
      for (int i = 0; i < nl[0]; i++){
        int s = k*nl[0]*nl[1] + j*nl[0] + i;
        tmp[s][0] = exponent * sin(_pi*x/length) * sin(_pi*y/length)* sin(_pi*z/length);
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
T accuracy(const std::vector<impl::StencilCell<T, M...>> &in1, std::vector<impl::StencilCell<T, M...>> &in2) {

  double err   = 0.0;
  double norm  = 0.0;

  for (int i = 0; i < in1.size(); i++){
    err  += (in1[i][0] - in2[i][0])*(in1[i][0] - in2[i][0]);
    norm += in1[i][0]*in1[i][0];
  }
  if constexpr (rel_error) return (T)sqrt(err/norm);
  else                     return (T)sqrt(err);
}
