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

//could be even more "generic", e.g. ND stencil
template <typename T, StencilType ST, int D = 3>
struct Generic3DStencil {

  using F = typename field_mapper<T, D>::type;

  F       out;
  const F in;

  // Stencil params here:
  const T c0;
  const T c1;

  Generic3DStencil(std::vector<T> &f1, const std::vector<T> &f2, const T c0_, const T c1_, const std::array<int, D> dims) :
	out(f1, dims),
	in (const_cast<std::vector<T>&>(f2), dims),
	c0(c0_),
	c1(c1_) { }

  template <int dir = 0>  
  inline T add_face_neighbors(const std::array<int, D> &x, const int i) {

    auto accum = in.template operator()<face_shifts[dir]> (x, i);
    
    if constexpr (dir < (2*D-1)) {
      constexpr int next_dir = dir + 1;
      return (accum + add_face_neighbors<next_dir>(x, i));  
    } else {//end recursion
      return accum;    
    }  
  } 
  
  inline T add_edge_neighbors(const std::array<int, D> &x, const int i) {
    return(in.template operator()<Shift::ShiftXp1, Shift::ShiftYp1> (x, i)+
           in.template operator()<Shift::ShiftXp1, Shift::ShiftYm1> (x, i)+
           in.template operator()<Shift::ShiftXm1, Shift::ShiftYp1> (x, i)+
           in.template operator()<Shift::ShiftXm1, Shift::ShiftYm1> (x, i)+
           in.template operator()<Shift::ShiftYp1, Shift::ShiftZp1> (x, i)+
           in.template operator()<Shift::ShiftYp1, Shift::ShiftZm1> (x, i)+
           in.template operator()<Shift::ShiftYm1, Shift::ShiftZp1> (x, i)+
           in.template operator()<Shift::ShiftYm1, Shift::ShiftZm1> (x, i)+
           in.template operator()<Shift::ShiftZp1, Shift::ShiftXp1> (x, i)+
           in.template operator()<Shift::ShiftZp1, Shift::ShiftXm1> (x, i)+
           in.template operator()<Shift::ShiftZm1, Shift::ShiftXp1> (x, i)+
           in.template operator()<Shift::ShiftZm1, Shift::ShiftXm1> (x, i));
  }

  inline T add_corner_neighbors(const std::array<int, D> &x, const int i) {
    return(in.template operator()<Shift::ShiftXp1, Shift::ShiftYp1, Shift::ShiftZp1> (x, i)+
           in.template operator()<Shift::ShiftXp1, Shift::ShiftYp1, Shift::ShiftZm1> (x, i)+
           in.template operator()<Shift::ShiftXp1, Shift::ShiftYm1, Shift::ShiftZp1> (x, i)+
           in.template operator()<Shift::ShiftXp1, Shift::ShiftYm1, Shift::ShiftZm1> (x, i)+
           in.template operator()<Shift::ShiftXm1, Shift::ShiftYp1, Shift::ShiftZp1> (x, i)+
           in.template operator()<Shift::ShiftXm1, Shift::ShiftYp1, Shift::ShiftZm1> (x, i)+
           in.template operator()<Shift::ShiftXm1, Shift::ShiftYm1, Shift::ShiftZp1> (x, i)+
           in.template operator()<Shift::ShiftXm1, Shift::ShiftYm1, Shift::ShiftZm1> (x, i));
  }


  typename std::enable_if<D == 3, void>::type operator()(const int i){//
     std::array<int, D> x{0};

     in.Indx2Coord(x, i);

     T res = c0*in[i]+c1*add_face_neighbors(x,i);

     if constexpr (ST == StencilType::FaceEdgeCentered || ST == StencilType::FaceEdgeCornerCentered) {
       res += c1*(0.5*add_edge_neighbors(x,i));
       //       
       if constexpr (ST == StencilType::FaceEdgeCornerCentered)  res += c1*(_1div3_*add_corner_neighbors(x,i));
     }

     out[i] = res;

     return;
  }
};
