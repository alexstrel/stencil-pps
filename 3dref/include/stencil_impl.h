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

//could be even more "generic", e.g. ND stencil
template <typename T, StencilType ST, int D = 3>
struct Generic3DStencil {

  using F = typename field_mapper<T, D>::type;

  F v;//data accessor

  // Stencil params here:
  const T c0;
  const T c1;
  const T c2;
  const T c3;  

  Generic3DStencil(std::vector<T> &latt, const T c0_, const T c1_, const std::array<int, D> dims) :
	v(latt.data(),dims),
	c0(c0_),
	c1(c1_),
	c2(c1_* 0.5),
	c3(c1_ * _1div3_){ }

  Generic3DStencil(const T c0_, const T c1_, const std::array<int, D> dims) :
	v(dims),
	c0(c0_),
	c1(c1_),
	c2(c1_* 0.5),
	c3(c1_ * _1div3_){ }


  void SetData(std::vector<T> &latt) {v.Set(latt);}
  
  template <int dir = 0>  
  inline T add_face_neighbors(const std::array<int, D> &x, const int i) {

    auto neigh = v.template operator()<shifts[dir]> (x, i);
    
    if constexpr (dir < (2*D-1)) {
      constexpr int next_dir = dir + 1;
      return (neigh + add_face_neighbors<next_dir>(x, i));  
    } else {//end recursion
      return neigh;    
    }  
  } 
  
  inline T add_edge_neighbors(const std::array<int, D> &x, const int i) {
    return(v.template operator()<Shift::ShiftXp1, Shift::ShiftYp1> (x, i)+
           v.template operator()<Shift::ShiftXp1, Shift::ShiftYm1> (x, i)+
           v.template operator()<Shift::ShiftXm1, Shift::ShiftYp1> (x, i)+
           v.template operator()<Shift::ShiftXm1, Shift::ShiftYm1> (x, i)+
           v.template operator()<Shift::ShiftYp1, Shift::ShiftZp1> (x, i)+
           v.template operator()<Shift::ShiftYp1, Shift::ShiftZm1> (x, i)+
           v.template operator()<Shift::ShiftYm1, Shift::ShiftZp1> (x, i)+
           v.template operator()<Shift::ShiftYm1, Shift::ShiftZm1> (x, i)+
           v.template operator()<Shift::ShiftZp1, Shift::ShiftXp1> (x, i)+
           v.template operator()<Shift::ShiftZp1, Shift::ShiftXm1> (x, i)+
           v.template operator()<Shift::ShiftZm1, Shift::ShiftXp1> (x, i)+
           v.template operator()<Shift::ShiftZm1, Shift::ShiftXm1> (x, i));
  }

  inline T add_corner_neighbors(const std::array<int, D> &x, const int i) {
    return(v.template operator()<Shift::ShiftXp1, Shift::ShiftYp1, Shift::ShiftZp1> (x, i)+
           v.template operator()<Shift::ShiftXp1, Shift::ShiftYp1, Shift::ShiftZm1> (x, i)+
           v.template operator()<Shift::ShiftXp1, Shift::ShiftYm1, Shift::ShiftZp1> (x, i)+
           v.template operator()<Shift::ShiftXp1, Shift::ShiftYm1, Shift::ShiftZm1> (x, i)+
           v.template operator()<Shift::ShiftXm1, Shift::ShiftYp1, Shift::ShiftZp1> (x, i)+
           v.template operator()<Shift::ShiftXm1, Shift::ShiftYp1, Shift::ShiftZm1> (x, i)+
           v.template operator()<Shift::ShiftXm1, Shift::ShiftYm1, Shift::ShiftZp1> (x, i)+
           v.template operator()<Shift::ShiftXm1, Shift::ShiftYm1, Shift::ShiftZm1> (x, i));
  }



  inline typename std::enable_if<D <= 4, T>::type operator()(const int i){
    std::array<int, D> x{0};
    
    v.Indx2Coord(x, i);    
    
    T res = c0*v[i] + c1*add_face_neighbors(x,i);

    if      constexpr (ST == StencilType::FaceEdgeCentered       && D > 1) res += c2*add_edge_neighbors(x,i);
    else if constexpr (ST == StencilType::FaceEdgeCornerCentered && D > 2) res += c2*add_edge_neighbors(x,i)+c3*add_corner_neighbors(x,i);

    return res;
  }
  
  typename std::enable_if<D <= 4, double>::type operator()(const T &in_v, T &out_v){//
  
     const int i = &in_v - &v[0];
     
     out_v = this->operator()(i);     

     return sqr(static_cast<double>(out_v - in_v));
  }  
  
};
