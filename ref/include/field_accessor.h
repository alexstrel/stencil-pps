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

    //works for both nvc++ and g++
    //template recursion (still ugly ...)
    template<Shift shift, Shift... other_shifts>
    inline constexpr T& GetNeighbor(const std::array<int,D> &x, const int i) const {
      int j;
      //currently unsafe: no check on dimensionality (should be done during stencil inst.)
      if     (shift == Shift::ShiftXp1) j = x[0] == N[0]-1 ? i : i+1;
      else if(shift == Shift::ShiftXm1) j = x[0] == 0      ? i : i-1;
      else if(shift == Shift::ShiftYp1) j = x[1] == N[1]-1 ? i : i+N[0];
      else if(shift == Shift::ShiftYm1) j = x[1] == 0      ? i : i-N[0];
      else if(shift == Shift::ShiftZp1) j = x[2] == N[2]-1 ? i : i+NxNy;
      else if(shift == Shift::ShiftZm1) j = x[2] == 0      ? i : i-NxNy;
      else j = i;

      if constexpr (sizeof...(other_shifts) != 0) {
        return GetNeighbor<other_shifts...>(x, j);
      } else {
        return v[j];//return the value and stop recursion..
      }
    }

    T& operator[](const int i) const { return v[i];}

    template<Shift face_shift, Shift... other_shifts>
    T& operator()(const std::array<int,D> &x, const int i) const {
      if constexpr ((sizeof...( other_shifts)  > 3)) {
        printf("Number of shifts is not supported.\n"); 
        exit(-1);
      }

      if constexpr (face_shift == Shift::NoShift) return v[i]; //no shifts, direct access
      else                                        return GetNeighbor<face_shift, other_shifts...>(x,i); //shifted access

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
