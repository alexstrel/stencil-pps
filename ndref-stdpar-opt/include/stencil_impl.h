#pragma once

#include <common.h>
#include <field.h>
#include <enums.h>
#include <field_accessor.h>

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

template <typename data_tp, int D, int blockx, int blocky, int blockz, ArithmeticTp... coeffs>
class GenericNDStencilArg {
  
  public:
    using T = data_tp::value_type;
    using S = FieldArgs<D>;
    using F = FieldAccessor<data_tp, S>;
    //
    static constexpr int Dims{D};
    static constexpr int Dims_ = std::max(D,3);
    //
    static constexpr int Bx{blockx};
    static constexpr int By{D > 1 ? blocky : 1};    
    static constexpr int Bz{D > 2 ? blockz : 1};    
    //
    S accessor_args;

    F out;//  
    F in;//stencil source , but not const!

    // Stencil params here:
    const std::array<int, Dims_> block_offsets;//only 3d blocking is supported
    const std::array<T, sizeof...(coeffs)> c;

    GenericNDStencilArg(std::vector<data_tp> &out_, const std::vector<data_tp> &in_,  const std::array<int, D> dims, const coeffs& ...c_) :
	accessor_args(dims),    
    	out(out_, accessor_args),
	in (const_cast<std::vector<data_tp>&>(in_), accessor_args),
	block_offsets([d=dims]()->decltype(auto) {
	        std::array<int, 3> b{blockx, blocky, blockz};
	        std::array<int, Dims_> tmp{d[0] / b[0], 1 ,1};//Nx/Bx, NxNy/(BxBy), NxNyNz/(BxByBz), NxNyNzNt/(BxByBz), etc..
	        //
                for(int i = 1; i < Dims_; i++) {
		  tmp[i] = tmp[i-1]*( i < 3 ? d[i] / b[i] : d[i]);
                } return tmp;} ()),	
	c{c_...} { 
    }   
    //
    void Swap() { out.swap(in); }
    
    inline decltype(auto) Indx2BlockCoord(const int &i) const {
       //
       std::array<int, Dims_> x;

       x[0] = i; // return it for 1D domain, otherwise use it also as temp.

       if constexpr (D > 2) {
         // First, compute higher dim coords:
#pragma unroll
         for (int j = D - 1; j > 1; j--) {
           x[j] = x[0] / block_offsets[j-1];
           x[0] = (x[0] - x[j]*block_offsets[j-1]);
         }
       }
       //
       if constexpr (D > 1) {
         x[1] = x[0] / block_offsets[0];
         x[0] = x[0] - x[1]*block_offsets[0];
       }

       return x;
     }    
};

//could be even more "generic", e.g. ND stencil
template <StencilTp ST, typename Args>
class GenericNDStencil {
  using Tp = typename Args::T;
  //
  static constexpr int D{Args::Dims};//real dimensions,
  static constexpr int D_{Args::Dims_};//needed for 3D blocking,
  // 
  static constexpr int BlckX{Args::Bx};
  static constexpr int BlckY{Args::By};
  static constexpr int BlckZ{Args::Bz};
  static constexpr int BlckV{Args::Bx*Args::By*Args::Bz};  

  const Args& arg;
  
  public:

  constexpr GenericNDStencil(const Args &arg) : arg(arg) {}	
  
  template < int init_dir = 0, int base_dir = 2, int dir = 4 >
  inline Tp add_corner_neighbors(const int face_type, const std::array<int, D_> &x, const int i, const int j) {
    //
    auto is_curr_dir_bndry =  check_stencil_bndry<init_dir, base_dir, dir>(face_type);
    //
    auto neigh  = is_curr_dir_bndry == 0 ? arg.in.template operator()<shifts[init_dir], shifts[base_dir], shifts[dir] > (i,j) :  arg.in.get_bndry_term<dir>(face_type,x,i,j);

    if constexpr (dir % 2 == 0) {
      return (neigh + add_corner_neighbors<init_dir, base_dir, dir+1>(face_type, x, i, j));
    } else if constexpr (base_dir % 2 == 0) {
      return (neigh + add_corner_neighbors<init_dir, base_dir+1>(face_type, x, i, j));
    } else if constexpr (init_dir % 2 == 0) {
      return (neigh + add_corner_neighbors<init_dir+1>(face_type, x, i, j));
    }
    // 
    return neigh;
  }

  template <int base_dir = 0, int dir = 2>
  inline Tp add_edge_neighbors(const int face_type, const std::array<int, D_> &x, const int i, const int j) {
    //
    auto is_curr_dir_bndry =  check_stencil_bndry<base_dir, dir>(face_type);
    //
    auto neigh  = is_curr_dir_bndry == 0 ? arg.in.template operator()<shifts[base_dir], shifts[dir] > (i,j) : arg.in.get_bndry_term<dir>(face_type,x,i,j);
    //
    if        constexpr ( dir % 2 == 0) {
      return (neigh + add_edge_neighbors<base_dir, dir+1>(face_type, x, i, j));
    } else if constexpr (base_dir < (2*D-1)) {//
      constexpr int next_dir      = 2*(((base_dir+1) / 2 + 1) % D);
      return (neigh + add_edge_neighbors<base_dir+1, next_dir>(face_type, x, i, j));
    }
    //end recursion
    return neigh;
  }

  template <int dir = 0>
  inline Tp add_neighbors(int &face_type, const std::array<int, D_> &x, const int i, const int j) {
    //
    int current_face_type = arg.in.check_face_type<dir>(x, j);
    //
    auto neigh  = current_face_type == 0 ? arg.in.template operator()<shifts[dir]> (i,j) : arg.in.get_bndry_term<dir>(current_face_type,x,i,j);
    // Update face type information:
    if constexpr (ST == StencilTp::FaceEdgeCentered or ST == StencilTp::FaceEdgeCornerCentered){
      face_type = face_type | current_face_type;
    }
    //
    if constexpr (dir < (2*D-1)) {
      return (neigh + add_neighbors<dir+1>(face_type, x, i, j));
    }
    // 
    return neigh;
  }

  inline typename std::enable_if<D <= 3, void>::type operator()(const int k){
    //
    const int j = k % arg.in.stencil_grid_size;     
    // Compute base coord:
    int  i = (k / arg.in.stencil_grid_size) * BlckV; 
    //
    auto x = arg.in.Indx2Coord(i);
    // Keep 3D base point
    const std::array<int, 3> y{x[0], x[1], x[2]}; 
#pragma unroll 
    for(int bz = 0; bz < BlckZ; bz++ ) {//perform prefetching if necessary
      x[2] = y[2] + bz;
#pragma unroll 
      for(int by = 0; by < BlckY; by++ ) {//perform prefetching if necessary
        x[1] = y[1] + by;
#pragma unroll 
        for(int bx = 0; bx < BlckX; bx++ ) {//perform prefetching if necessary
          //
          x[0] = y[0] + bx;
          i    = arg.in.Coord2Indx(x);
          int face_type = 0;
          //
          auto res = arg.c[0]*arg.in.template operator()<Shift::NoShift> (i, j) + arg.c[1]*add_neighbors(face_type,x,i,j);

          if      constexpr (ST == StencilTp::FaceEdgeCentered  && D > 1)      res += arg.c[2]*add_edge_neighbors(face_type,x,i,j);
          else if constexpr (ST == StencilTp::FaceEdgeCornerCentered && D > 2) res += arg.c[2]*add_edge_neighbors(face_type,x,i,j)+arg.c[3]*add_corner_neighbors(face_type,x,i,j);

          arg.out[i][j] = res;          
        }
      }      
    }
    
    return;
  }

  typename std::enable_if<D <= 4, double>::type operator()(Tp &out_v, const Tp &in_v){//
  
     auto sqr = [](double x) { return x * x; };
  
     const int i = &in_v - &arg.in[0];
     
     this->operator()(i);          

     return sqr(static_cast<double>(out_v - in_v));
  }

};  
