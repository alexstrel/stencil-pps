#pragma once

#include <enums.h>
#include <field_accessor.h>

constexpr double _1div3_ = 1.0 / 3.0;

template<typename T>
void debug_7pt_stencil
(std::vector<T> &out, const std::vector<T> &in, const int Lx, const int Ly, const int Lz, const float C0, const float C1)
{
#pragma omp parallel num_threads(4)
  for (int k = 0; k < Lz; k++){
#pragma omp for schedule(dynamic)
    for (int j = 0; j < Ly; j++){
#pragma omp simd
      for (int i = 0; i < Lx; i++) {
        const int s = i + j*Lx + k*Lx*Ly;
        const T xp1 = i == Lx-1  ? 0.0 : in[s + 1];
        const T xm1 = i == 0     ? 0.0 : in[s - 1];
        const T yp1 = j == Ly-1  ? 0.0 : in[s + Lx];
        const T ym1 = j == 0     ? 0.0 : in[s - Lx];
        const T zp1 = k == Lz-1  ? 0.0 : in[s + Lx*Ly];
        const T zm1 = k == 0     ? 0.0 : in[s - Lx*Ly];

        out[s] = C0*in[s]+C1*(xp1+xm1+yp1+ym1+zp1+zm1);
      }
    }
  }
  return;
}

auto sqr = [](double x) { return x * x; };

#if 1
template <KernelType type> inline bool doBulk() {
  switch (type) {
    case Exterior_Kernel_X:
    case Exterior_Kernel_Y:
    case Exterior_Kernel_Z: return false;
    case Interior_Kernel:   return true;
  }

  return false;
}

template <KernelType type> inline bool doHalo(int dim = -1)
{
  switch (type) {
    case Exterior_Kernel_X: return dim == 0 || dim == -1 ? true : false;
    case Exterior_Kernel_Y: return dim == 1 || dim == -1 ? true : false;
    case Exterior_Kernel_Z: return dim == 2 || dim == -1 ? true : false;
    case Interior_Kernel:   return false;
  }
  return false;
}

#endif

//could be even more "generic", e.g. ND stencil
template <typename T, StencilType ST, int D = 3>
struct Generic3DStencil {

  using F = typename field_mapper<T, D>::type;

  F v;//data accessor

  // Stencil params here:
  const std::array<T, 4> c;

  Generic3DStencil(std::vector<T> &latt, const T c0_, const T c1_, const std::array<int, D> dims) :
	v(latt.data(),dims),
	c{c0_, c1_, 0.5, 0.5*_1div3_}
        { }

  Generic3DStencil(const T c0_, const T c1_, const std::array<int, D> dims) :
	v(dims),
	c{c0_, c1_, 0.5, 0.5*_1div3_}
        { }

  void SetData(std::vector<T> &latt) {v.Set(latt);}

  template < int init_dir, int base_dir, int dir >
  inline T add_corner_neighbors(const std::array<int, D> &x, const int i) {
    //
    bool is_bndry = v.check_bndry<dir>(x);
    //
    auto neigh  = is_bndry == false ? v.template operator()<shifts[init_dir], shifts[base_dir], shifts[dir] > (x, i) :  v.get_bndry_term<dir>(x,i);
    // 
    if constexpr (dir < (2*D-1)) {//Note that each dim has two dirs
      constexpr int next_dir = dir + 1;
      return (neigh + add_corner_neighbors<init_dir, base_dir, next_dir>(x, i));
    }
    // 
    return neigh;
  }

  template <int base_dir, int dir >
  inline T add_edge_neighbors(const std::array<int, D> &x, const int i) {
    //
    bool is_bndry = v.check_bndry<dir>(x);
    //
    auto neigh  = is_bndry == false ? v.template operator()<shifts[base_dir], shifts[dir] > (x, i) : v.get_bndry_term<dir>(x,i);
    //
    if constexpr (ST == StencilType::FaceEdgeCornerCentered && (dir == 2 || dir == 3)) {
      if(is_bndry == false) {
	constexpr int next_dir = 4;//Always start with Z dirs
        neigh = neigh + c[3]*add_corner_neighbors<base_dir, dir, next_dir>(x,i);
      }
    }
    //
    if constexpr (dir % 2 == 0) {//Note that each dim has two dirs
      constexpr int next_dir = dir + 1;
      return (neigh + add_edge_neighbors<base_dir, next_dir>(x, i));
    } 
    //end recursion
    return neigh;
  }

  template <int dir = 0>
  inline T add_neighbors(const std::array<int, D> &x, const int i) {
    //
    bool is_bndry = v.check_bndry<dir>(x);
    //
    auto neigh  = is_bndry == false ? v.template operator()<shifts[dir]> (x, i) : v.get_bndry_term<dir>(x,i);
    //
    if constexpr (ST == StencilType::FaceEdgeCentered || ST == StencilType::FaceEdgeCornerCentered) {	    
      if(is_bndry == false) { 
	constexpr int next_dir = 2*((dir / 2 + 1) % D);   
        neigh = neigh + c[2]*add_edge_neighbors<dir, next_dir>(x,i);
      }	
    }
    //
    if constexpr (dir < (2*D-1)) {//Note that each dim has two dirs
      constexpr int next_dir = dir + 1;
      return (neigh + add_neighbors<next_dir>(x, i));
    }
    // 
    return neigh;
  }

  inline typename std::enable_if<D <= 3, T>::type operator()(const int i){
    std::array<int, D> x{0};

    v.Indx2Coord(x, i);
    //
    auto res = c[0]*v.template operator()<Shift::NoShift> (x,i) + c[1]*add_neighbors(x,i);

    return res;
  }

  typename std::enable_if<D <= 4, double>::type operator()(T &out_v, const T &in_v){//
  
     const int i = &in_v - &v[0];
     
     out_v = this->operator()(i);     

     return sqr(static_cast<double>(out_v - in_v));
  }

};  
