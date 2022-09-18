#pragma once
#include <common.h>

const float _pi  = acosf(-1.f);
const float _2pi = (2.f*_pi);

template<typename Policy,  typename T, typename Allocator>
void create_field(const Policy &p,
                  std::vector<T, Allocator> &out,
                  const std::array<int, 3> nl,
                  const T kappa,
                  const T length,
                  const T time) {
                  
  std::vector<T> tmp(out.size());                  
  
  const std::array<T, 3> dnl{length / (nl[0]+static_cast<T>(1.0)), length /(nl[1]+static_cast<T>(1.0)), length / (nl[2]+static_cast<T>(1.0))} ; 	
  
  const T exponent  = expf(-3.0f*kappa*_pi*_pi*time/(length*length));
  
  T z = dnl[2];
  
  for (int k = 0; k < nl[2]; k++){
    T y = dnl[1];
    for (int j = 0; j < nl[1]; j++){
      T x = dnl[0];
      for (int i = 0; i < nl[0]; i++){
        int s = k*nl[0]*nl[1] + j*nl[0] + i;
        tmp[s] = exponent * sinf(_pi*x/length) * sinf(_pi*y/length)* sinf(_pi*z/length);
        x += dnl[0];
      }
      y += dnl[1];
    }
    z += dnl[2];
  }
  
  std::copy(p, tmp.begin(), tmp.end(), out.begin());
  
  return;
}

template<typename T, typename Allocator, bool rel_error = true>
T accuracy(const std::vector<T, Allocator> &in1, std::vector<T, Allocator> &in2) {

  double err   = 0.0;
  double norm  = 0.0;

  for (int i = 0; i < in1.size(); i++){
    err  += (in1[i] - in2[i])*(in1[i] - in2[i]);
    norm += in1[i]*in1[i];
  }
  if constexpr (rel_error) return (T)sqrt(err/norm);
  else                     return (T)sqrt(err);
}
