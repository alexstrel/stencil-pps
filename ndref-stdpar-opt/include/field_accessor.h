#pragma once

#include <common.h>
#include <field.h>

template <typename StencilCell, int D=3>
class FieldAccessor{
  public:
    using T = typename StencilCell::value_type; 
    
    static constexpr int stencil_cell_size = StencilCell::cell_size;

  private :
    StencilCell *v;//no allocation

    const std::array<int, D> N;//domain size 

    const int NxNy;

  public :
    FieldAccessor(std::vector<StencilCell> &latt, const std::array<int, D> &dims) :
        v(latt.data()), 
        N(dims),
        NxNy(D > 1 ? dims[0]*dims[1] : 0) {
        }
        
    void Set(StencilCell* in) {v = in;}  
    StencilCell* Get() {return v;} 
    
    void swap(FieldAccessor &f){
      //
      StencilCell* tmp = f.Get();
      f.Set(v);
      v = tmp;
      
      return;
    }

    inline int GetDim(int dim) {return N[dim];}

    inline decltype(auto) GetDims() const {return N;} 

    template <Shift shift>
    inline int get_shift_dir() const {

      if constexpr        (shift == Shift::ShiftXp1) {
         return 0;
      } else if constexpr (shift == Shift::ShiftXm1) {
         return 1;
      } else if constexpr (shift == Shift::ShiftYp1) {
         return 2;
      } else if constexpr (shift == Shift::ShiftYm1) {
         return 3;
      } else if constexpr (shift == Shift::ShiftZp1) {
         return 4;
      } else if constexpr (shift == Shift::ShiftZm1) {
         return 5;
      }

      return -1;
    }	    

    template<int dir>
    inline constexpr int check_face_type(const std::array<int,D> &x, const int i) const {

      int domain_face_idx = 0;
      
      if        constexpr (dir == 0) {
        if (x[0] == N[0] - 1) domain_face_idx = 1; 
      } else if constexpr (dir == 1) {
        if (x[0] == 0       ) domain_face_idx = 2;
      } else if constexpr (dir == 2) {
        if (x[1] == N[1] - 1) domain_face_idx = 4;
      } else if constexpr (dir == 3) {
        if (x[1] == 0       ) domain_face_idx = 8;
      } else if constexpr (dir == 4) {
        if (x[2] == N[2] - 1) domain_face_idx =16;
      } else if constexpr (dir == 5) {
        if (x[2] == 0) domain_face_idx =32;
      }
      
      if(domain_face_idx != 0) {
        std::array<int, StencilCell::D> y{0};

        StencilCell::Indx2Coord(y, i);
        
        if constexpr (dir == 0) {
          if (domain_face_idx == 1 && y[0] != StencilCell::m[0] - 1) domain_face_idx = 0; 
        } else if constexpr (dir == 1) {
          if (domain_face_idx == 2 && y[0] != 0                    ) domain_face_idx = 0;
        } else if constexpr (dir == 2) {
          if (domain_face_idx == 4 && y[1] != StencilCell::m[1] - 1) domain_face_idx = 0;
        } else if constexpr (dir == 3) {
          if (domain_face_idx == 8 && y[1] != 0                    ) domain_face_idx = 0;
        } else if constexpr (dir == 4) {
         if (domain_face_idx == 16 && y[2] != StencilCell::m[2] - 1) domain_face_idx = 0;
        } else if constexpr (dir == 5) {
         if (domain_face_idx == 32 && y[2] != 0                    ) domain_face_idx = 0;
        }        
      }

      return domain_face_idx;
    }

    template<int dir>
    inline T get_bndry_term(const std::array<int, D> &x, const int j, const int i) const {
      return static_cast<T>(0.0);//Trivial BC
    }

    //works for both nvc++ and g++
    //template recursion (still ugly ...)
    template<Shift shift, Shift... other_shifts>
    inline constexpr int GetNeighborIdx(int j) const {
      //currently unsafe: no check on dimensionality (should be done during stencil inst.)
      if constexpr        (shift == Shift::ShiftXp1) {
         j += 1;
      } else if constexpr (shift == Shift::ShiftXm1) {
         j -= 1;
      } else if constexpr (shift == Shift::ShiftYp1) {
         j += N[0];
      } else if constexpr (shift == Shift::ShiftYm1) {
         j -= N[0];
      } else if constexpr (shift == Shift::ShiftZp1) {
         j += NxNy;
      } else if constexpr (shift == Shift::ShiftZm1) {
         j -= NxNy;
      }

      if constexpr (sizeof...(other_shifts) != 0) {
        return GetNeighborIdx<other_shifts...>(j);
      } 

      return j;//stop recursion and return the index of the neighbor
    }

    StencilCell& operator[](const int j) const { return v[j];}

    template<Shift face_shift, Shift... other_shifts>
    T operator()(const int j, const int i) const {
      //	      
      if constexpr ((sizeof...( other_shifts)  > 3)) {
        printf("Number of shifts is not supported.\n"); 
        exit(-1);
      }
      //
      if constexpr (face_shift == Shift::NoShift) return v[j][i]; //no shifts, direct access

      const auto k = GetNeighborIdx<face_shift, other_shifts...>(j); //shifted access

      return v[k][i];
    }

    //3d version
    inline void Indx2Coord(std::array<int, D> &x, const int &i) const {
      x[2] = i / NxNy;
      const int tmp = (i - x[2]*NxNy);
      x[1] = tmp / N[0];
      x[0] = tmp - x[1]*N[0];
    }
};
