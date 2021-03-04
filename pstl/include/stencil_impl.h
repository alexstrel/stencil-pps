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
template <typename T, int M, StencilType StTp, int D>
struct GenericNDStencil {

  using F  = typename field_mapper<T, M, D>::type; //type refers to std::array<T,M>
  using TM = std::array<T,M>;
  using ID = std::array<int,D>;

  F       out;
  const F in;

  // Stencil params here:
  const T c0;
  const T c1;

  GenericNDStencil(std::vector<TM> &f1, const std::vector<TM> &f2, const T c0_, const T c1_, const ID dims) :
	out(f1, dims),
	in (const_cast<std::vector<TM>&>(f2), dims),
	c0(c0_),
	c1(c1_){ }

  template<typename Args>
  GenericNDStencil(MDLattice<T, M, D, Args> &l) :
          out(l.V(), l.Extents()),
          in(l.Tmp1(), l.Extents()),
          c0(l.GetArgs().c0),
          c1(l.GetArgs().c1){ }

  inline T add_face_neighbors(const ID &x, const int j) {
    return(in.template operator()<Shift::ShiftXp1> (x, j)+
           in.template operator()<Shift::ShiftXm1> (x, j)+
           in.template operator()<Shift::ShiftYp1> (x, j)+
           in.template operator()<Shift::ShiftYm1> (x, j)+
           in.template operator()<Shift::ShiftZp1> (x, j)+
           in.template operator()<Shift::ShiftZm1> (x, j));
  }

  inline T add_edge_neighbors(const ID &x, const int j) {
    return(in.template operator()<Shift::ShiftXp1, Shift::ShiftYp1> (x, j)+
           in.template operator()<Shift::ShiftXp1, Shift::ShiftYm1> (x, j)+
           in.template operator()<Shift::ShiftXm1, Shift::ShiftYp1> (x, j)+
           in.template operator()<Shift::ShiftXm1, Shift::ShiftYm1> (x, j)+
           in.template operator()<Shift::ShiftYp1, Shift::ShiftZp1> (x, j)+
           in.template operator()<Shift::ShiftYp1, Shift::ShiftZm1> (x, j)+
           in.template operator()<Shift::ShiftYm1, Shift::ShiftZp1> (x, j)+
           in.template operator()<Shift::ShiftYm1, Shift::ShiftZm1> (x, j)+
           in.template operator()<Shift::ShiftZp1, Shift::ShiftXp1> (x, j)+
           in.template operator()<Shift::ShiftZp1, Shift::ShiftXm1> (x, j)+
           in.template operator()<Shift::ShiftZm1, Shift::ShiftXp1> (x, j)+
           in.template operator()<Shift::ShiftZm1, Shift::ShiftXm1> (x, j));
  }

  inline T add_corner_neighbors(const ID &x, const int j) {
    return(in.template operator()<Shift::ShiftXp1, Shift::ShiftYp1, Shift::ShiftZp1> (x, j)+
           in.template operator()<Shift::ShiftXp1, Shift::ShiftYp1, Shift::ShiftZm1> (x, j)+
           in.template operator()<Shift::ShiftXp1, Shift::ShiftYm1, Shift::ShiftZp1> (x, j)+
           in.template operator()<Shift::ShiftXp1, Shift::ShiftYm1, Shift::ShiftZm1> (x, j)+
           in.template operator()<Shift::ShiftXm1, Shift::ShiftYp1, Shift::ShiftZp1> (x, j)+
           in.template operator()<Shift::ShiftXm1, Shift::ShiftYp1, Shift::ShiftZm1> (x, j)+
           in.template operator()<Shift::ShiftXm1, Shift::ShiftYm1, Shift::ShiftZp1> (x, j)+
           in.template operator()<Shift::ShiftXm1, Shift::ShiftYm1, Shift::ShiftZm1> (x, j));
  }


  inline T compute_site_stencil(const ID &x, const int j) {
    //if constexpr (D != 3) { std::cout << "Stencil dimensionality is not supported." << std::endl; exit(-1); }

    switch(StTp) {
      case StencilType::FaceCentered  :
        return c0*in[j]+c1*add_face_neighbors(x,j);
        break;
      case StencilType::FaceEdgeCentered :
        return c0*in[j]+c1*(add_face_neighbors(x,j)+0.5*add_edge_neighbors(x,j));
        break;
      case StencilType::FaceEdgeCornerCentered :
        return c0*in[j]+c1*(add_face_neighbors(x,j)+0.5*add_edge_neighbors(x,j)+_1div3_*add_corner_neighbors(x,j));
        break;
      default :
        exit(-1);
        break;
    }
  }

  template <StencilPolicy spolicy = StencilPolicy::DefaultPolicy>
  typename std::enable_if<D == 3, void>::type operator()(const int j){//
     ID x{0};

     switch(spolicy) {
       case StencilPolicy::DefaultPolicy :
         {
           //compute base coordinate:
           int l = j*M;
           in.Indx2Coord(x, l);
#pragma unroll
           for(int i = 0; i < M; i++) {
             out[l] = compute_site_stencil(x, l);
	     l += 1; x[0] += 1;
           }
         }
         break;
       default :
         exit(-1);
         break;
     }

     return;
  }
};
