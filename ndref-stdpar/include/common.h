#pragma once

#include <concepts> 
#include <type_traits>

template <typename T>
concept ArithmeticTp = std::is_arithmetic_v<T>; 


