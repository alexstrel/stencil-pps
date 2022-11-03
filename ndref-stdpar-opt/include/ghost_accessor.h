#pragma once

#include <common.h>
#include <field.h>

template<int... M>
consteval decltype(auto) compute_cell_offsets() {

  std::array<int, sizeof... (M)> m{M...};
  //
  std::array<int, sizeof... (M)> offsets{M...};
  //
  int prev_offset = 1; 
  //
  for (int i = 0; i < m.size(); i++) { //NB: does not check for Mx=1, My=1, Mz=1
    offsets[i]  = prev_offset*m[i];//m[0], m[0]*m[1], m[0]*m[1]*m[2]
    prev_offset = offsets[i];
  }

  return offsets;
}

template <int D, int... M>
class GhostArgs{
  public:
    //
    static constexpr std::array<int, sizeof...(M)> m{M...};
    static constexpr std::array<int, sizeof...(M)> cell_offsets{compute_cell_offsets<M...>()};    
    static constexpr int cell_dim{sizeof...(M)};    

    const std::pair<std::array<int, D>, std::array<int, D>> params;

    GhostArgs(const std::array<int, D> &dims) :
	params([&d=dims]()->std::pair<std::array<int, D>, std::array<int, D>> {
                std::array<int, D> nm1{d[0]-1};		
                std::array<int, D> offsets{d[0]};//Nx, NxNy, NxNyNz, ... etc
                std::array<int, D> strides{d[0]-1};//Nx-1, NxNy-Nx, NxNyNz-NxNy, ... etc
		//
                for(int i = 1; i < D; i++) {
                  nm1[i]     = d[i] - 1;
                  offsets[i] = offsets[i-1]*d[i];
                  strides[i] = offsets[i] - offsets[i-1];
                } return std::make_pair(nm1, strides);} ()) {}

    inline decltype(auto) get_rangem1() const { return params.first;}
    inline decltype(auto) get_strides() const { return params.second;}
    
};

template <typename FAccessor, typename Arg>
class GhostAccessor{
  public:
    using T = typename FAccessor::T; //?

    static constexpr int D                 = FAccessor::D;
    static constexpr int stencil_cell_size = FAccessor::stencil_cell_size;    

  private :
    const FAccessor& v;    
    const Arg& args;

  public :
    GhostAccessor(const FAccessor &u, const Arg &args) : 
        v(u),
        args(args){ 
        }

    inline static decltype(auto) CellIndx2CellCoord(const int &i) {
//    inline decltype(auto) CellIndx2CellCoord(const int &i) const {        
      //
      constexpr auto& cell_offsets = Arg::cell_offsets;
      constexpr auto& cell_dim     = Arg::cell_dim;
      //      
      std::array<int, cell_dim> y{i};          
#if 0 // potential compiler bug here!
#pragma unroll         
      for (int j = (cell_dim-1); j >= 1; j--) {                  
        y[j] = y[0] / cell_offsets[j-1];
        y[0] = (y[0] - y[j]*cell_offsets[j-1]);
      }
#else // this works perfect!
      y[2] = y[0] / cell_offsets[1];
      y[0] = y[0] - y[2]*cell_offsets[1];
      // 
      y[1] = y[0] / cell_offsets[0];
      y[0] = y[0] - y[1] * cell_offsets[0];      
#endif
       
      return y;     
    }    

    template<int dir>
    inline constexpr int check_face_type(const std::array<int,D> &x, const int i) const {

      constexpr int num_dirs = 2*D;
      constexpr int shift    = num_dirs;

      const auto& Nm1 = args.get_rangem1();

      int domain_face_idx    = 0;

      if        constexpr (dir == 0) {
        if (x[0] == Nm1[0]  ) domain_face_idx = 1; 
      } else if constexpr (dir == 1) {
        if (x[0] == 0       ) domain_face_idx = 2;
      } else if constexpr (dir == 2) {
        if (x[1] == Nm1[1]  ) domain_face_idx = 4;
      } else if constexpr (dir == 3) {
        if (x[1] == 0       ) domain_face_idx = 8;
      } else if constexpr (dir == 4) {
        if (x[2] == Nm1[2]  ) domain_face_idx =16;
      } else if constexpr (dir == 5) {
        if (x[2] == 0       ) domain_face_idx =32;
      }
      
      if constexpr (stencil_cell_size > 1) {
        if(domain_face_idx != 0) {

          const auto y = CellIndx2CellCoord(i);
        
          if constexpr (dir == 0) {
            if (domain_face_idx ==  1 && y[0] < Arg::m[0] - 1) domain_face_idx <<= shift; 
          } else if constexpr (dir == 1) {
            if (domain_face_idx ==  2 && y[0] > 0            ) domain_face_idx <<= shift;
          } else if constexpr (dir == 2) {
            if (domain_face_idx ==  4 && y[1] < Arg::m[1] - 1) domain_face_idx <<= shift;
          } else if constexpr (dir == 3) {
            if (domain_face_idx ==  8 && y[1] > 0            ) domain_face_idx <<= shift;
          } else if constexpr (dir == 4) {
            if (domain_face_idx == 16 && y[2] < Arg::m[2] - 1) domain_face_idx <<= shift;
          } else if constexpr (dir == 5) {
            if (domain_face_idx == 32 && y[2] > 0            ) domain_face_idx <<= shift;
          }        
        }
      }
      return domain_face_idx;
    }

    template<Shift shift, Shift... other_shifts>
    inline T get_bndry_term(const int face_type, const std::array<int, D> &x, const int j, const int i) const {
      //
      const auto& strides = args.get_strides();	    

      if constexpr (stencil_cell_size > 1){
        if constexpr        (shift == Shift::ShiftXp1) {
          if (face_type & 64  ) {
            const int k = j-strides[0];// Nm1[0];
            const int l = i+1;
            if constexpr (sizeof...(other_shifts) != 0) {
	      return get_bndry_term<other_shifts...>(face_type, x, k, l); 
	    }
	    return v.template operator()<Shift::NoShift>(k, l);
          }
        } else if constexpr (shift == Shift::ShiftXm1) {
          if (face_type & 128 ) {
            const int k = j+strides[0];
            const int l = i-1;
	    if constexpr (sizeof...(other_shifts) != 0) {
              return get_bndry_term<other_shifts...>(face_type, x, k, l);
            }
	    return v.template operator()<Shift::NoShift>(k, l);
          }
        } else if constexpr (shift == Shift::ShiftYp1) {
          if (face_type & 256 ) {
            const int k = j-strides[1];//NxNymNx;
            const int l = i+Arg::m[0];
            if constexpr (sizeof...(other_shifts) != 0) {
              return get_bndry_term<other_shifts...>(face_type, x, k, l);
            }
	    return v.template operator()<Shift::NoShift>(k, l);
          }
        } else if constexpr (shift == Shift::ShiftYm1) {
          if (face_type & 512 ){
            const int k = j+strides[1];//NxNymNx;
            const int l = i-Arg::m[0];
            if constexpr (sizeof...(other_shifts) != 0) {
              return get_bndry_term<other_shifts...>(face_type, x, k, l);
            }
	    return v.template operator()<Shift::NoShift>(k, l);
          }
        } else if constexpr (shift == Shift::ShiftZp1) {
          if (face_type & 1024){
            const int k = j-strides[2];//NxNyNzmNxNy;
            const int l = i+Arg::m[0]*Arg::m[1];
            if constexpr (sizeof...(other_shifts) != 0) {
              return get_bndry_term<other_shifts...>(face_type, x, k, l);
            }
	    return v.template operator()<Shift::NoShift>(k, l);
          }
        } else if constexpr (shift == Shift::ShiftZm1) {
          if (face_type & 2048){
            const int k = j+strides[2];//NxNyNzmNxNy;
            const int l = i-Arg::m[0]*Arg::m[1];
            if constexpr (sizeof...(other_shifts) != 0) {
              return get_bndry_term<other_shifts...>(face_type, x, k, l);
            }
  	    return v.template operator()<Shift::NoShift>(k, l);
          }
        }  
      }      
      // For "true" boundary:  
      return static_cast<T>(0.0);//Trivial BC
      
    }

};
