#pragma once

#include <concepts> 
#include <type_traits>

// Simple arithmetic type
template <typename T>
concept ArithmeticTp = requires{
  std::is_arithmetic_v<T>;
}

// Simple complex type
template <typename T>
concept ComplexTp    = requires (T t) {
  requires ArithmeticTp<decltype(t.real())>;
  requires ArithmeticTp<decltype(t.imag())>;  
}


