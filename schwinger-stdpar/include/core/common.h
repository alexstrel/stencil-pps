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

#include <limits>

constexpr int min_int                = std::numeric_limits<int>::min();
constexpr std::size_t min_sizet      = std::numeric_limits<std::size_t>::min();

//some defs:
constexpr std::size_t invalid_spin  = min_sizet;
constexpr std::size_t invalid_color = min_sizet;
constexpr std::size_t invalid_dir   = min_sizet;

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

// Simple generic arithmetic type
template <typename T>
concept ArithmeticTp = FloatTp<T> or ComplexTp<T>;


// Generic container type:
template <typename T, typename = std::void_t<>> class is_container_type : public std::false_type { };

template <typename T> class is_container_type< T, std::void_t<decltype(std::declval<T>().begin()), 
                                                              decltype(std::declval<T>().end()), 
                                                              decltype(std::declval<T>().data()), 
                                                              decltype(std::declval<T>().size())> > : public std::true_type { };
                                                              
template <typename T> constexpr bool is_container_type_v = is_container_type<T>::value;                                                              

template <typename T>
concept GenericContainerTp = is_container_type_v<T>; 

// More specific container type (allocator aware containers, e.g., std::vector etc.):
//
template <typename T, typename = std::void_t<>> class is_allocated_type : public std::false_type { };

template <typename T> class is_allocated_type< T, std::void_t<typename T::allocator_type>> : public std::true_type { };

template <typename T> constexpr bool is_allocated_type_v = is_allocated_type<T>::value;

// Dynamic container type:
template <typename T>
concept ContainerTp = GenericContainerTp<T> and is_allocated_type_v<T>;

template <typename T>
concept PMRContainerTp = ContainerTp<T> and std::is_same< typename T::allocator_type,  std::pmr::polymorphic_allocator<typename T::value_type> >::value;

template <typename T>
concept ContainerViewTp = GenericContainerTp<T> and (not is_allocated_type_v<T>);

// Iterator type
template <typename T>
concept IteratorTp = std::random_access_iterator<T>;


