#pragma once

#include <common.h>
#include <field.h>

template <int D, int... M>
class FieldArgs{
  public:
    static constexpr int Dims{D};
    static constexpr int cell_size{impl::get_cell_size<M...>()};    
    //
    const std::array<int, D> offsets;

    FieldArgs(const std::array<int, D> &dims) :
	offsets([&d=dims]()->std::array<int, D>{
                std::array<int, D> nm1{d[0]-1};		
                std::array<int, D> offsets{d[0]};//Nx, NxNy, NxNyNz, ... etc
		//
                for(int i = 1; i < D; i++) {
                  nm1[i]     = d[i] - 1;
                  offsets[i] = offsets[i-1]*d[i];
                } return offsets;} ()) {}
};

template <typename StencilCell, typename Arg>
class FieldAccessor{
  public:
    using T = typename StencilCell::value_type; 
    static constexpr int D{Arg::Dims};
    static constexpr int E{std::max(Arg::Dims, 3)};    
    static constexpr int stencil_cell_size{Arg::cell_size};    

  private :
    StencilCell *v;//no allocation
    
    const Arg& args;

  public :
    FieldAccessor(std::vector<StencilCell> &latt, const Arg &args) :
        v(latt.data()), 
        args(args){ 
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

    //works for both nvc++ and g++
    //template recursion (still ugly ...)
    template<Shift shift, Shift... other_shifts>
    inline constexpr int GetNeighborIdx(int j) const {
 
      const auto& offsets = args.offsets;    

      //currently unsafe: no check on dimensionality (should be done during stencil inst.)
      if constexpr        (shift == Shift::ShiftXp1) {
         j += 1;
      } else if constexpr (shift == Shift::ShiftXm1) {
         j -= 1;
      } else if constexpr (shift == Shift::ShiftYp1) {
         j += offsets[0];//Nx
      } else if constexpr (shift == Shift::ShiftYm1) {
         j -= offsets[0];//Nx
      } else if constexpr (shift == Shift::ShiftZp1) {
         j += offsets[1];//NxNy
      } else if constexpr (shift == Shift::ShiftZm1) {
         j -= offsets[1];//NxNy
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

    inline decltype(auto) Indx2Coord(const int &i) const {
       //
       const auto& offsets = args.offsets;
       //
       std::array<int, E> x{i};// return it for 1D domain, otherwise use it also as temp.

       if constexpr (D > 1) {
#pragma unroll
         for (int j = D - 1; j >= 1; j--) {
           x[j] = x[0] / offsets[j-1];
           x[0] = (x[0] - x[j]*offsets[j-1]);
         }
       }

       return x;
     }
     
     inline int Coord2Indx(const std::array<int, E> &x) const {
       //
       const auto& offsets = args.get_offsets();       
       //
       int i = x[0];
#pragma unroll
       for (int j = 1; j < D; j++) {
         i += x[j]*offsets[j-1];
       }

       return i;
     }
};
