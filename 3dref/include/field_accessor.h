#pragma once

template <typename T, int D=3>
class FieldAccessor{
  private :
    T *v;//no allocation

    const std::array<int, D> N; 

    const int NxNy;

  public :
    FieldAccessor(std::vector<T> &latt, const std::array<int, D> &dims) :
        v(latt.data()), 
	N(dims),
        NxNy(D > 1 ? dims[0]*dims[1] : 0) {
        }

    FieldAccessor(const std::array<int, D> &dims) :
        v(nullptr), 
	N(dims),
        NxNy(D > 1 ? dims[0]*dims[1] : 0) {
        }
        
    void Set(std::vector<T> &latt) {v = latt.data();}    

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

    template <Shift shift1, Shift shift2>
    inline int get_shift_dir() const {
      return get_shift_dir<shift2>();	    
    }

    template <Shift shift1, Shift shift2, Shift shift3>
    inline int get_shift_dir() const {
      return get_shift_dir<shift3>();
    }

    template<int dir>
    inline bool check_bndry(const std::array<int,D> &x) const {

      constexpr int bndry_dir = dir % 2;
      constexpr int dim       = dir / 2;

      int bndry_coord = 0;
      //
      if constexpr (bndry_dir == 0) {
        bndry_coord = N[dim] - 1;
      }
      //
      return (x[dim] == bndry_coord ? true : false);
    }
    

    template<int dir>
    inline T get_bndry_term(const std::array<int, D> &x, const int i) const {
      return static_cast<T>(0.0);//Trivial BC
    }

    //works for both nvc++ and g++
    //template recursion (still ugly ...)
    template<Shift shift, Shift... other_shifts>
    inline constexpr int GetNeighborIdx(const std::array<int,D> &x, int i) const {
      //currently unsafe: no check on dimensionality (should be done during stencil inst.)
      if constexpr        (shift == Shift::ShiftXp1) {
         i += 1;
      } else if constexpr (shift == Shift::ShiftXm1) {
         i -= 1;
      } else if constexpr (shift == Shift::ShiftYp1) {
         i += N[0];
      } else if constexpr (shift == Shift::ShiftYm1) {
         i -= N[0];
      } else if constexpr (shift == Shift::ShiftZp1) {
         i += NxNy;
      } else if constexpr (shift == Shift::ShiftZm1) {
         i -= NxNy;
      }

      if constexpr (sizeof...(other_shifts) != 0) {
        return GetNeighborIdx<other_shifts...>(x, i);
      } 

      return i;//stop recursion and return the index of the neighbor
    }

    T& operator[](const int i) const { return v[i];}

    template<Shift face_shift, Shift... other_shifts>
    T operator()(const std::array<int,D> &x, const int i) const {
      //	      
      if constexpr ((sizeof...( other_shifts)  > 3)) {
        printf("Number of shifts is not supported.\n"); 
        exit(-1);
      }
      //
      if constexpr (face_shift == Shift::NoShift) return v[i]; //no shifts, direct access

      const auto j = GetNeighborIdx<face_shift, other_shifts...>(x,i); //shifted access

      return v[j];
    }

    //3d version
    inline void Indx2Coord(std::array<int, D> &x, const int &i) const {
      x[2] = i / NxNy;
      const int tmp = (i - x[2]*NxNy);
      x[1] = tmp / N[0];
      x[0] = tmp - x[1]*N[0];
    }
};

template <typename T, int D> struct field_mapper {
};

template <int D> struct field_mapper<double, D> {
  typedef FieldAccessor<double, D> type;
};

template <int D> struct field_mapper<float, D> {
  typedef FieldAccessor<float, D> type;
};
