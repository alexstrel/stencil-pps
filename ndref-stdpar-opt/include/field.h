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

template <ArithmeticTp T, int... M>
class StencilCell {
   public :
     using value_type = T;
     //local dimensions
     static constexpr int D{sizeof...(M)};
     static constexpr std::array<int, D> m{M...};
     static constexpr int cell_size{getsize<M...>()};

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
     inline static void Indx2Coord(std::array<int, D> &x, const int &i) {
       if constexpr (D == 3) {
         constexpr int NxNy = m[0]*m[1];
         x[2] = i / NxNy;
         const int tmp = (i - x[2]*NxNy);
         x[1] = tmp / m[0];
         x[0] = tmp - x[1]*m[0];
       }
     }
   
};

}//end namespace impl
