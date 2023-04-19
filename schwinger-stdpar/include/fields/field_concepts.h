#pragma once

#include <core/common.h>

#include <memory_resource>

// Generic Spinor Field type:
template <typename T>
concept GenericSpinorFieldTp = requires{
  requires GenericContainerTp<typename T::container_tp>;
  requires (T::nSpin  == 1ul or T::nSpin  == 2ul or T::nSpin  == 4ul);
  requires (T::nColor == 1ul or T::nColor == 3ul);
  requires (T::nDir   == invalid_dir); 
};

// Generic Gauge Field type:
template <typename T>
concept GenericGaugeFieldTp = requires{
  requires GenericContainerTp<typename T::container_tp>;
  requires (T::nSpin  == invalid_spin);
  requires (T::nColor == 1ul or T::nColor == 3ul);
  requires (T::nDir   >= 2ul and T::nDir  <= 4ul);
};

// Generic Field type :
template <typename T>
concept GenericFieldTp    = GenericSpinorFieldTp<T> or GenericGaugeFieldTp<T>;

// Allocated field type
template <typename T>
concept FieldTp           = GenericFieldTp<T> and is_allocated_type_v<typename T::container_tp>;

// Reference field type
template <typename T>
concept FieldViewTp       = GenericFieldTp<T> and (not is_allocated_type_v<typename T::container_tp>);

// Allocated field type
template <typename T>
concept SpinorFieldTp     = GenericSpinorFieldTp<T> and is_allocated_type_v<typename T::container_tp>;

// Reference field type
template <typename T>
concept SpinorFieldViewTp = GenericSpinorFieldTp<T> and (not is_allocated_type_v<typename T::container_tp>);

// PMR spinor field type
template <typename T>
concept PMRSpinorFieldTp = SpinorFieldTp<T> and (std::is_same_v<typename T::container_tp::allocator_type, std::pmr::polymorphic_allocator< typename T::container_tp::value_type > >);

// Allocated field type
template <typename T>
concept GaugeFieldTp      = GenericGaugeFieldTp<T> and is_allocated_type_v<typename T::container_tp>;

template <typename T>
concept GaugeFieldViewTp  = GenericGaugeFieldTp<T> and (not is_allocated_type_v<typename T::container_tp>);

// Spinor block Field concepts
template <typename T>
concept BlockSpinorFieldTp     = SpinorFieldTp<T> or (ContainerTp<T> and SpinorFieldTp< typename std::remove_pointer< decltype( std::declval<T>().data() ) >::type >); 

template <typename T>
concept BlockSpinorFieldViewTp = SpinorFieldViewTp<T> or (ContainerViewTp<T> and SpinorFieldViewTp< typename std::remove_pointer< decltype( std::declval<T>().data() ) >::type >); 

template <typename T>
concept PMRBlockSpinorFieldTp  = ContainerTp<T> and PMRSpinorFieldTp< typename std::remove_pointer< decltype( std::declval<T>().data() ) >::type >;


