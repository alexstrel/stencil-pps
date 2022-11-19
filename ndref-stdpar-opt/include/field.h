#pragma once

#include <common.h>

namespace impl 
{

template<int... M>
consteval int get_cell_size() {

  std::array<int, sizeof... (M)> r{M...};
  int s = 1;

  for (int &i : r)  s = s*i;

  return s;
}

template <ArithmeticTp T, int... M>
class StencilCell {
   public :
     using value_type = T;

     T data[get_cell_size<M...>()];

     StencilCell() = default;
     StencilCell(const StencilCell<T, M...> &) = default;
     StencilCell(StencilCell<T, M...> &&)      = default;
   
     //basic accessors
     constexpr T &operator[](const int i) { return data[i]; }
     constexpr const T &operator[](const int i) const { return data[i]; }

     constexpr int size() const { return get_cell_size<M...>(); }   

     auto operator=(const StencilCell&) -> StencilCell& = default;
     auto operator=(StencilCell&&     ) -> StencilCell& = default;
};

}//end namespace impl
