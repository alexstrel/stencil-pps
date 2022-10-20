#pragma once

#include <concepts> 
#include <type_traits>

// Simple arithmetic type
template <typename T>
concept FloatTp = requires{
  std::is_floating_point_v<T>;
};

// Simple complex type
template <typename T>
concept ComplexTp    = requires (T t) {
  requires FloatTp<decltype(t.real())>;
  requires FloatTp<decltype(t.imag())>;  
};


// Simple genetic arithmetic type
template <typename T>
concept ArithmeticTp = FloatTp<T> || ComplexTp<T>;
