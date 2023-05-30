#pragma once

#include <limits>

constexpr int min_int      = std::numeric_limits<int>::min();
constexpr std::size_t min_sizet      = std::numeric_limits<std::size_t>::min();

//some defs:
constexpr std::size_t invalid_spin  = min_sizet;
constexpr std::size_t invalid_color = min_sizet;
constexpr std::size_t invalid_dir   = min_sizet;

//using a-la QUDA terminology 
enum class FieldSiteSubset {
  FullSiteSubset         = 2,  
  ParitySiteSubset       = 1,  
  InvalidSiteSubset      = min_int  
};

enum class FieldParity {
  EvenFieldParity        =  0, 
  OddFieldParity         =  1, 
  InvalidFieldParity     =  min_int  
};

enum class FieldType {
  ScalarFieldType    =  1,
  VectorFieldType    =  2,
  SpinorFieldType    =  3,
  InvalidFieldType   =  min_int
};

enum class FieldOrder {
  LexFieldOrder     =  1,
  EOFieldOrder      =  2,
  InvalidFieldOrder =  min_int
};


