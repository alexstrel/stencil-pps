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

template<int dir, int... other_dirs>
inline constexpr int check_stencil_bndry(const int face_idx, int face_type = 0) {

  int is_bndry = 0;

  if        constexpr (dir == 0) {
    is_bndry = face_idx & 1;
  } else if constexpr (dir == 1) {
    is_bndry = face_idx & 2;
  } else if constexpr (dir == 2) {
    is_bndry = face_idx & 4;
  } else if constexpr (dir == 3) {
    is_bndry = face_idx & 8;
  } else if constexpr (dir == 4) {
    is_bndry = face_idx & 16;
  } else if constexpr (dir == 5) {
    is_bndry = face_idx & 32;
  }

  face_type = face_type | is_bndry;

  if constexpr (sizeof...(other_dirs) != 0) {
    return check_stencil_bndry<other_dirs...>(face_idx, face_type);
  }

  return face_type;
}


//could be even more "generic", e.g. ND stencil
template <typename T, StencilType ST, int D = 3>
struct Generic3DStencil {

  using F = typename field_mapper<T, D>::type;

  F v;//data accessor

  // Stencil params here:
  const std::array<T, 4> c;

  Generic3DStencil(std::vector<T> &latt, const T c0_, const T c1_, const std::array<int, D> dims) :
	v(latt.data(),dims),
	c{c0_, c1_, c1_* 0.5, c1_*_1div3_}
        { }

  Generic3DStencil(const T c0_, const T c1_, const std::array<int, D> dims) :
	v(dims),
	c{c0_, c1_, c1_* 0.5, c1_*_1div3_}
        { }

  void SetData(std::vector<T> &latt) {v.Set(latt);}

  template < int init_dir = 0, int base_dir = 2, int dir = 4 >
  inline T add_corner_neighbors(const int face_type, const std::array<int, D> &x, const int i) {
    //
    auto is_curr_dir_bndry =  check_stencil_bndry<init_dir, base_dir, dir>(face_type);
    //
    auto neigh  = is_curr_dir_bndry == 0 ? v.template operator()<shifts[init_dir], shifts[base_dir], shifts[dir] > (i) :  v.get_bndry_term<dir>(x,i);

    if constexpr (dir % 2 == 0) {
      return (neigh + add_corner_neighbors<init_dir, base_dir, dir+1>(face_type, x, i));
    } else if constexpr (base_dir % 2 == 0) {
      return (neigh + add_corner_neighbors<init_dir, base_dir+1>(face_type, x, i));
    } else if constexpr (init_dir % 2 == 0) {
      return (neigh + add_corner_neighbors<init_dir+1>(face_type, x, i));
    }
    // 
    return neigh;
  }

  template <int base_dir = 0, int dir = 2>
  inline T add_edge_neighbors(const int face_type, const std::array<int, D> &x, const int i) {
    //
    auto is_curr_dir_bndry =  check_stencil_bndry<base_dir, dir>(face_type);
    //
    auto neigh  = is_curr_dir_bndry == 0 ? v.template operator()<shifts[base_dir], shifts[dir] > (i) : v.get_bndry_term<dir>(x,i);
    //
    if        constexpr ( dir % 2 == 0) {
      return (neigh + add_edge_neighbors<base_dir, dir+1>(face_type, x, i));
    } else if constexpr (base_dir < (2*D-1)) {//
      constexpr int next_dir      = 2*(((base_dir+1) / 2 + 1) % D);
      return (neigh + add_edge_neighbors<base_dir+1, next_dir>(face_type, x, i));
    }
    //end recursion
    return neigh;
  }

  template <int dir = 0>
  inline T add_neighbors(int &face_type, const std::array<int, D> &x, const int i) {
    //
    auto current_face_type = v.check_face_type<dir>(x);
    //
    auto neigh  = current_face_type == 0 ? v.template operator()<shifts[dir]> (i) : v.get_bndry_term<dir>(x,i);
    // Update face type information:
    if constexpr (ST == StencilType::FaceEdgeCentered or ST == StencilType::FaceEdgeCornerCentered){
      face_type = face_type | current_face_type;
    }
    //
    if constexpr (dir < (2*D-1)) {
      return (neigh + add_neighbors<dir+1>(face_type, x, i));
    }
    // 
    return neigh;
  }

  inline typename std::enable_if<D <= 3, T>::type operator()(const int i){
    std::array<int, D> x{0};

    v.Indx2Coord(x, i);
    //
    int face_type = 0;
    //
    auto res = c[0]*v.template operator()<Shift::NoShift> (i) + c[1]*add_neighbors(face_type, x, i);

    if      constexpr (ST == StencilType::FaceEdgeCentered  && D > 1)      res += c[2]*add_edge_neighbors(face_type, x,i);
    else if constexpr (ST == StencilType::FaceEdgeCornerCentered && D > 2) res += c[2]*add_edge_neighbors(face_type, x,i)+c[3]*add_corner_neighbors(face_type, x,i);

    return res;
  }

  typename std::enable_if<D <= 4, double>::type operator()(T &out_v, const T &in_v){//
  
     const int i = &in_v - &v[0];
     
     out_v = this->operator()(i);     

     return sqr(static_cast<double>(out_v - in_v));
  }

};  
