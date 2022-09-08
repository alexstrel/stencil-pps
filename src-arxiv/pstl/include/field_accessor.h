#pragma once

template <typename T, int M,  int D>
class FieldAccessor{
  private :
    T *v;//raw pointer

    const std::array<int, D> N; //N[0]=>Nx, N[1]=>Ny, N[2]=>Nz etc
//we perform communications in upto 4 dim
    const int NxNy;
    const int NxNymNx;
    const int NxNyNz;
    const int NxNyNzmNxNy;
    const int NxNyNzNt;
    const int NxNyNzNtmNxNyNz;

  public :
    FieldAccessor(std::vector<std::array<T, M>> &latt, const std::array<int, D> &dims) :
        v(reinterpret_cast<T*>(latt.data())), 
	N(dims),
        NxNy(D > 1 ? dims[0]*dims[1] : 0),
        NxNymNx(D > 1 ? dims[0]*dims[1]-dims[0] : 0),
        NxNyNz(D > 2 ? dims[0]*dims[1]*dims[2] : 0),
        NxNyNzmNxNy(D > 2 ? NxNyNz-NxNy : 0),
        NxNyNzNt(D > 3 ? dims[0]*dims[1]*dims[2]*dims[3] : 0),
        NxNyNzNtmNxNyNz(D > 3 ? NxNyNzNt-NxNyNz : 0){  }

    //works for both nvc++ and g++
    //template recursion (still ugly ...)
    template<Shift shift, Shift... other_shifts>
    inline constexpr int GetNeighborIdx(const std::array<int,D> &x, int i) const {
      //currently unsafe: no check on dimensionality (should be done during stencil inst.)
      if constexpr        (shift == Shift::ShiftXp1) {
         if (x[0] == (N[0]-1))  return -1;
         i += 1;
      } else if constexpr (shift == Shift::ShiftXm1) {
         if (x[0] == 0       )  return -1;
         i -= 1;
      } else if constexpr (shift == Shift::ShiftYp1) {
         if (x[1] == (N[1]-1))  return -1;
         i += N[0];
      } else if constexpr (shift == Shift::ShiftYm1) {
         if (x[1] == 0       )  return -1;
         i -= N[0];
      } else if constexpr (shift == Shift::ShiftZp1) {
	 if (x[2] == (N[2]-1))  return -1;
         i += NxNy;
      } else if constexpr (shift == Shift::ShiftZm1) {
         if (x[2] == 0       )  return -1;
         i -= NxNy;
      } else if constexpr (shift == Shift::ShiftTp1) {
	 if (x[3] == (N[3]-1))  return -1;
         i += NxNyNz;
      } else if constexpr (shift == Shift::ShiftTm1) {
         if (x[3] == 0       )  return -1;
         i -= NxNyNz;
      }

      if constexpr (sizeof...(other_shifts) != 0) {
        return GetNeighborIdx<other_shifts...>(x, i);
      } else {
        return i;//stop recursion and return an index of the neighbor
      }
    }


    T& operator[](const int j) const { return v[j];}

    template<Shift face_shift, Shift... other_shifts>
    T operator()(const std::array<int,D> &x, const int i) const {
      if constexpr (((sizeof...( other_shifts) != 0) && D == 1) ||
                    ((sizeof...( other_shifts)  > 1) && D == 2) ||
                    ((sizeof...( other_shifts)  > 2) && D == 3) ||
                    ((sizeof...( other_shifts)  > 3) && D  > 4) ) {
        printf("Number of shifts is not supported.\n"); exit(-1);
      }
      
      if constexpr (face_shift == Shift::NoShift) return v[i]; //no shifts, direct access
      const auto j = GetNeighborIdx<face_shift, other_shifts...>(x,i); //shifted access

      return (j == -1 ? static_cast<T>(0.0) : v[j]);      
    }
    //3d version
    inline void Indx2Coord(std::array<int, D> &x, const int &i) const {
      x[2] = i / NxNy;
      const int tmp = (i - x[2]*NxNy);
      x[1] = tmp / N[0];
      x[0] = tmp - x[1]*N[0];
    }
};

template <typename T, int M, int D> struct field_mapper {
};

template <int M, int D> struct field_mapper<double, M, D> {
  typedef FieldAccessor<double, M, D> type;
};

template <int M, int D> struct field_mapper<float, M, D> {
  typedef FieldAccessor<float, M, D> type;
};
