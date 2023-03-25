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

// Simple genetic arithmetic type
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

// More specific container type (memory-owner containers, i.e., excl. container adapters):
//
template <typename T, typename = std::void_t<>> class is_allocated_type : public std::false_type { };

template <typename T> class is_allocated_type< T, std::void_t<typename T::allocator_type> > : public std::true_type { };

template <typename T> constexpr bool is_allocated_type_v = is_allocated_type<T>::value;

template <typename T>
concept AllocatedContainerTp = is_container_type_v<T> and is_allocated_type_v<T>;

template <typename T>
concept ReferenceContainerTp = is_container_type_v<T> and (not is_allocated_type_v<T>);

// Iterator type
template <typename T>
concept IteratorTp = std::random_access_iterator<T>;

// Generic Spinor Field type:
template <typename T>
concept GenericSpinorFieldTp = requires{
  requires GenericContainerTp<typename T::container_tp>;
  requires (T::nSpin  >= 1ul);
  requires (T::nColor >= 1ul);
  requires (T::nDir   == invalid_dir); 
};

// Generic Gauge Field type:
template <typename T>
concept GenericGaugeFieldTp = requires{
  requires GenericContainerTp<typename T::container_tp>;
  requires (T::nSpin  == invalid_spin);
  requires (T::nColor >= 1ul);
  requires (T::nDir   >= 2ul);
};

// Generic Field type :
template <typename T>
concept GenericFieldTp = GenericSpinorFieldTp<T> or GenericGaugeFieldTp<T>;

// Allocated field type
template <typename T>
concept AllocatedFieldTp = GenericFieldTp<T> and is_allocated_type_v<typename T::container_tp>;

// Reference field type
template <typename T>
concept ReferenceFieldTp = GenericFieldTp<T> and (not is_allocated_type_v<typename T::container_tp>);

// Block Field
template <typename T>
concept ReferenceBlockFieldTp = (ReferenceContainerTp<T> and ReferenceFieldTp<typename T::value_tp>); 

// Allocated field type
template <typename T>
concept AllocatedSpinorFieldTp = GenericSpinorFieldTp<T> and is_allocated_type_v<typename T::container_tp>;

// Allocated field type
template <typename T>
concept AllocatedGaugeFieldTp = GenericGaugeFieldTp<T> and is_allocated_type_v<typename T::container_tp>;

// Reference field type
template <typename T>
concept ReferenceSpinorFieldTp = GenericSpinorFieldTp<T> and (not is_allocated_type_v<typename T::container_tp>);

template <typename T>
concept ReferenceGaugeFieldTp = GenericGaugeFieldTp<T> and (not is_allocated_type_v<typename T::container_tp>);


