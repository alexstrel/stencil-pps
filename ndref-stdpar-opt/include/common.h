#pragma once

#include <concepts> 
#include <type_traits>

template <typename T>
concept ArithmeticTp = requires { 
  std::is_arithmetic_v<T>; 
};


