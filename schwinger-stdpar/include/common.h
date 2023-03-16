#pragma once

#include <concepts> 
#include <type_traits>
#include <span>

#include <ranges>
#include <vector>
#include <iostream>
#include <complex>
#include <algorithm>

#include <random>

#include <experimental/mdspan>

namespace stdex = std::experimental;

// Simple arithmetic type
template <typename T>
concept FloatTp = std::is_floating_point_v<T>;

// Simple complex type
template <typename T>
concept ComplexTp    = requires {
  requires std::is_floating_point_v<decltype(std::declval<T>().real())>;
  requires std::is_floating_point_v<decltype(std::declval<T>().imag())>;  
};

// Simple genetic arithmetic type
template <typename T>
concept ArithmeticTp = FloatTp<T> || ComplexTp<T>;


// Generic container type:
template <typename T, typename = std::void_t<>> class is_container_type : public std::false_type { };

template <typename T> class is_container_type< T, std::void_t<decltype(std::declval<T>().begin()), 
                                                              decltype(std::declval<T>().end()), 
                                                              decltype(std::declval<T>().data()), 
                                                              decltype(std::declval<T>().size())> > : public std::true_type { };
                                                              
template <typename T> constexpr bool is_container_type_v = is_container_type<T>::value;                                                              

template <typename T>
concept GenericContainerTp = is_container_type_v<T>; 


