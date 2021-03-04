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
    inline constexpr T& GetNeighbor(const std::array<int,D> &x, const int i) const {
      int j;
      //currently unsafe: no check on dimensionality (should be done during stencil inst.)
      if     (shift == Shift::ShiftXp1) j = x[0] == N[0]-1 ? i : i+1;
      else if(shift == Shift::ShiftXm1) j = x[0] == 0      ? i : i-1;
      else if(shift == Shift::ShiftYp1) j = x[1] == N[1]-1 ? i : i+N[0];
      else if(shift == Shift::ShiftYm1) j = x[1] == 0      ? i : i-N[0];
      else if(shift == Shift::ShiftZp1) j = x[2] == N[2]-1 ? i : i+NxNy;
      else if(shift == Shift::ShiftZm1) j = x[2] == 0      ? i : i-NxNy;
      else if(shift == Shift::ShiftTp1) j = x[3] == N[3]-1 ? i : i+NxNyNz;
      else if(shift == Shift::ShiftTm1) j = x[3] == 0      ? i : i-NxNyNz;      
      else j = i;

      if constexpr (sizeof...(other_shifts) != 0) {//dpcpp !
        return GetNeighbor<other_shifts...>(x, j);
      } else {
        return v[j];//return the value and stop recursion..
      }
    }

    T& operator[](const int j) const { return v[j];}

    template<Shift face_shift, Shift... other_shifts>
    T& operator()(const std::array<int,D> &x, const int i) const {
      if constexpr (((sizeof...( other_shifts) != 0) && D > 1) ||
                    ((sizeof...( other_shifts)  > 1) && D > 2) ||
                    ((sizeof...( other_shifts)  > 2) && D > 3) ||
                    ((sizeof...( other_shifts)  > 3) && D > 4) ) {
        std::cout << "Number of shifts is not supported." << std::endl; exit(-1);
      }

      if(face_shift == Shift::NoShift) return v[i]; //no shifts, direct access
      else                             return GetNeighbor<face_shift, other_shifts...>(x,i); //shifted access
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
