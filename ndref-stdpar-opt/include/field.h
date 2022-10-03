#pragma once

#include <common.h>

namespace impl 
{

template<int... M>
consteval int getsize() {

  std::array<int, sizeof... (M)> r{M...};
  int s = 1;

  for (int &i : r)  s = s*i;

  return s;
}

template<int... M>
consteval decltype(auto) compute_global_offsets() {

  std::array<int, sizeof... (M)> m{M...};
  //
  std::array<int, sizeof... (M)> offsets{M...};
  //
  int i = 0;
  //
#if 0  
  offsets[0] = m[0]; 
  //
  for (int &m_ : m) { 
    offests[i+1] = offsets[i]*m_;//!
    i += 1;
  }
#else
  int prev_offset = 1; 
  //
  for (int &m_ : m) { 
    offests[i]  = prev_offset*m_;//m[0], m[0]*m[1], m[0]*m[1]*m[2]
    prev_offset = offests[i];
    i += 1;
  }
#endif
  return offsets;
}	

template <ArithmeticTp T, int... M>
class StencilCell {
   public :
     using value_type = T;

     static constexpr std::array<int, sizeof...(M)> m{M...};
     static constexpr int cell_size{getsize<M...>()};
     static constexpr std::array<int, sizeof...(M)> offsets{compute_global_offsets<M...>()};     

     T data[cell_size];

     StencilCell() = default;
     StencilCell(const StencilCell<T, M...> &) = default;
     StencilCell(StencilCell<T, M...> &&)      = default;
   
     //basic accessors
     constexpr T &operator[](const int i) { return data[i]; }
     constexpr const T &operator[](const int i) const { return data[i]; }

     constexpr int size() const { return getsize<M...>(); }   

     auto operator=(const StencilCell&) -> StencilCell& = default;
     auto operator=(StencilCell&&     ) -> StencilCell& = default;
     
   //private:
     inline static decltype(auto) Indx2Coord(const int &i) {
       //
       std::array<int, sizeof...(M)> x;     
     
       x[0] = i; // return for 1D domain, otherwise use also as temp.
           
       if constexpr (sizeof...(M) > 2) {
         // First, compute higher dim coords:
#pragma unroll         
         for (int j = (sizeof...(M)-1); j > 1; j--) {
           x[j] = x[0] / offsets[j-1];
           x[0] = (x[0] - x[j]*offsets[j-1]);
         }
       }
       //
       if constexpr (sizeof...(M) > 1) {
         x[1] = x[0] / offsets[0];
         x[0] = x[0] - x[1]*offsets[0];       
       } 
       
       return x;     
     }
   
};

}//end namespace impl
