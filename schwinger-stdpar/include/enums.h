#pragma once

#include <limits>

constexpr int min_int      = std::numeric_limits<int>::min();

//some defs:
constexpr int invalid_spin  = min_int;
constexpr int invalid_color = min_int;
constexpr int invalid_dir   = min_int;

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

