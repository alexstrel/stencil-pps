#pragma once
#include <common.h> //for concepts

template<ArithmeticTp T, int d = 3>
class NDField {
  private :

  std::vector<T> v;

  public:
    NDField(const int vol) : v{vol} {
    }

    ~NDField() {}//?
    
    auto& V() { return v;}

    template<typename Policy>
    void Set(const Policy &p, const T &&val){
      std::fill(policy, v.begin(), v.end(), val);
    }

    template<typename Policy>
    void Set(const Policy &p, const T &&val){
      std::fill(p, v.begin(), v.end(), val);
    }

    template<typename Policy>
    void Set(const Policy &p,                   
             const std::array<int, 3> nl,
             const T kappa,
             const T length,
             const T time){
             
      std::vector<Float> tmp(v.size());
      
      const std::array<T, 3> dnl{length / (nl[0]+static_cast<T>(1.0)), length /(nl[1]+static_cast<T>(1.0)), length / (nl[2]+static_cast<T>(1.0))} ; 	
      const T exponent  = exp(-3.0*kappa*_pi*_pi*time/(length*length));
      T z = dnl[2];
      for (int k = 0; k < nl[2]; k++){
        T y = dnl[1];
        for (int j = 0; j < nl[1]; j++){
          T x = dnl[0];
          for (int i = 0; i < nl[0]; i++){
            int s = k*nl[0]*nl[1] + j*nl[0] + i;
            out[s] = exponent * sin(_pi*x/length) * sin(_pi*y/length)* sin(_pi*z/length);
            x += dnl[0];
          }
          y += dnl[1];
        }
        z += dnl[2];
      }
      std::copy(p, tmp.begin(), tmp.end(), v.begin());
    }
};
