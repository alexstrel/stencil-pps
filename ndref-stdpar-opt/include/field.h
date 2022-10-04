#pragma once

#include <common.h>

namespace impl 
{

template<int... M>
consteval int get_grid_size() {

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
  int prev_offset = 1; 
  //
  for (int &m_ : m) { 
    offests[i]  = prev_offset*m_;//m[0], m[0]*m[1], m[0]*m[1]*m[2]
    prev_offset = offests[i];
    i += 1;
  }

  return offsets;
}	

template <ArithmeticTp T, int... M>
class StencilGrid {
   public :
     using value_type = T;

     static constexpr std::array<int, sizeof...(M)> m{M...};
     static constexpr int grid_size{get_grid_size<M...>()};
     static constexpr std::array<int, sizeof...(M)> offsets{compute_global_offsets<M...>()};     

     T data[grid_size];

     StencilGrid() = default;
     StencilGrid(const StencilGrid<T, M...> &) = default;
     StencilGrid(StencilGrid<T, M...> &&)      = default;
   
     //basic accessors
     constexpr T &operator[](const int i) { return data[i]; }
     constexpr const T &operator[](const int i) const { return data[i]; }

     constexpr int size() const { return get_grid_size<M...>(); }   

     auto operator=(const StencilGrid&) -> StencilGrid& = default;
     auto operator=(StencilGrid&&     ) -> StencilGrid& = default;
     
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
