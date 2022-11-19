#pragma once

#include <common.h>
#include <field.h>
#include <enums.h>
#include <field_accessor.h>
#include <ghost_accessor.h>


template<int dir, int... other_dirs>
inline constexpr int check_stencil_bndry(const int face_idx, int face_type = 0) {

  int is_bndry = 0;

  if        constexpr (dir == 0) {
    is_bndry = face_idx & 65;// 1 or 64 (inner cell)
  } else if constexpr (dir == 1) {
    is_bndry = face_idx & 130;//2 or 128 
  } else if constexpr (dir == 2) {
    is_bndry = face_idx & 260;//4 or 256
  } else if constexpr (dir == 3) {
    is_bndry = face_idx & 520;//8 or 512
  } else if constexpr (dir == 4) {
    is_bndry = face_idx & 1040;//16 or 1024
  } else if constexpr (dir == 5) {
    is_bndry = face_idx & 2080;//32 or 2048
  }

  face_type = face_type | is_bndry;

  if constexpr (sizeof...(other_dirs) != 0) {
    return check_stencil_bndry<other_dirs...>(face_idx, face_type);
  }

  return face_type;
}

<<<<<<< HEAD
template <typename data_tp, int D, int blockx, int blocky, int blockz, ArithmeticTp... coeffs>
=======
template <typename data_tp, int D, int blockx, int blocky, int blockz, int... M>
>>>>>>> develop
class GenericNDStencilArg {
  
  public:
    using T = data_tp::value_type;
    
    using U = FieldArgs<D, M...>;
    using F = FieldAccessor<data_tp, U>;
    //
    using V = GhostArgs<D, M...>;
    using G = GhostAccessor<F, V>;    
    //
    //
    static constexpr int Dims{D};
    static constexpr int EDims = std::max(D,3);
    //
    static constexpr int Bx{blockx};
    static constexpr int By{D > 1 ? blocky : 1};    
    static constexpr int Bz{D > 2 ? blockz : 1};    
    //
    U accessor_args;

    F out;//  
    F in;//stencil source , but not const!
    //
    V ghost_args;
    
    G in_ghost;    

    // Stencil params here:
    const std::array<T, 4> c;

    GenericNDStencilArg(std::vector<data_tp> &out_, const std::vector<data_tp> &in_,  const std::array<int, D> dims, const std::array<T, 4>& c_) :
	accessor_args(dims),    
    	out(out_, accessor_args),
	in (const_cast<std::vector<data_tp>&>(in_), accessor_args),
	ghost_args(dims),
	in_ghost (in, ghost_args),	    	
	c{c_} { 
    }   
    //
    void Swap() { out.swap(in); }//swap data pointers only!
};

//could be even more "generic", e.g. ND stencil
template <StencilTp ST, typename Args>
class GenericNDStencil {
  using Tp = typename Args::T;
  //
  static constexpr int D{Args::Dims};//real dimensions,
  static constexpr int E{Args::EDims};//needed for 3D blocking,
  // 
  static constexpr int BlockX{Args::Bx};
  static constexpr int BlockY{Args::By};
  static constexpr int BlockZ{Args::Bz};
  static constexpr int BlockV{Args::Bx*Args::By*Args::Bz};  

  const Args& arg;
  
  public:

  constexpr GenericNDStencil(const Args &arg) : arg(arg) {}	
  
  template < int first_dir = 0, int second_dir = 2, int third_dir = 4 >
  inline Tp add_corner_neighbors(const int face_type, const std::array<int, E> &x, const int i, const int j) {
    //
    auto is_curr_dir_bndry =  check_stencil_bndry<first_dir, second_dir, third_dir>(face_type);
    //
    auto neigh  = is_curr_dir_bndry == 0 ? arg.in.template operator()<shifts[first_dir], shifts[second_dir], shifts[third_dir] > (i,j) :  arg.in_ghost.template get_bndry_term<true, shifts[first_dir], shifts[second_dir], shifts[third_dir]>(face_type,x,i,j);

    if constexpr (third_dir % 2 == 0) {
      return (neigh + add_corner_neighbors<first_dir, second_dir, third_dir+1>(face_type, x, i, j));
    } else if constexpr (second_dir % 2 == 0) {
      return (neigh + add_corner_neighbors<first_dir, second_dir+1>(face_type, x, i, j));
    } else if constexpr (first_dir % 2 == 0) {
      return (neigh + add_corner_neighbors<first_dir+1>(face_type, x, i, j));
    }
    // 
    return neigh;
  }

  template <int first_dir = 0, int second_dir = 2>
  inline Tp add_edge_neighbors(const int face_type, const std::array<int, E> &x, const int i, const int j) {
    //
    auto is_curr_dir_bndry =  check_stencil_bndry<first_dir, second_dir>(face_type);
    //
    auto neigh  = is_curr_dir_bndry == 0 ? arg.in.template operator()<shifts[first_dir], shifts[second_dir] > (i,j) : arg.in_ghost.template get_bndry_term<true, shifts[first_dir], shifts[second_dir]>(face_type,x,i,j);
    //
    if        constexpr ( second_dir % 2 == 0) {
      return (neigh + add_edge_neighbors<first_dir, second_dir+1>(face_type, x, i, j));
    } else if constexpr (first_dir < (2*D-1)) {//
      constexpr int next_second_dir      = 2*(((first_dir+1) / 2 + 1) % D);
      return (neigh + add_edge_neighbors<first_dir+1, next_second_dir>(face_type, x, i, j));
    }
    //end recursion
    return neigh;
  }

  template <int dir = 0>
  inline Tp add_neighbors(int &face_type, const std::array<int, E> &x, const int i, const int j) {
    //
    auto current_face_type = arg.in_ghost.check_face_type<dir>(x, j);
    //
    auto neigh  = current_face_type == 0 ? arg.in.template operator()<shifts[dir]> (i,j) : arg.in_ghost.template get_bndry_term<true, shifts[dir]>(current_face_type,x,i,j);
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

  inline typename std::enable_if<D <= 3, void>::type operator()(const int l){
    //
    const int j = l % arg.in.stencil_cell_size;     
    // Compute base coord:
    const int i0= (l / arg.in.stencil_cell_size) * BlockV; 
    //
    auto x = arg.in.Indx2Coord(i0);
    // Keep 3D base point
    const std::array<int, 3> y{x[0], x[1], x[2]}; 
#pragma unroll 
    for(int bz = 0; bz < BlockZ; bz++ ) {//perform prefetching if necessary
      x[2] = y[2] + bz;
#pragma unroll 
      for(int by = 0; by < BlockY; by++ ) {//perform prefetching if necessary
        x[1] = y[1] + by;
#pragma unroll 
        for(int bx = 0; bx < BlockX; bx++ ) {//perform prefetching if necessary
          //
          x[0] = y[0] + bx;
          
          const auto i  = arg.in.Coord2Indx(x);

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
