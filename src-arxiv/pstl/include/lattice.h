#pragma once

#include <enums.h>

template<int D>
class MDLatticeParam {
  private :
    static constexpr int Ndim = D;
    std::array<int, D> nd;

    BCType bc;
    OrderType order;

    mutable bool parity;//false for the full lattice

    int vol;
    int volCB;

  public:
    explicit MDLatticeParam() : nd{0}, bc{BCType::Undefined}, order{OrderType::Undefined}, parity(false), vol(0), volCB(0){}

    MDLatticeParam(const std::array<int, D> nd_, BCType bc_ = BCType::Dirichlet, OrderType order_ = OrderType::Lexicographic, bool parity = false) : nd{nd_}, bc{bc_}, order(order_), parity(parity) {
      vol = 1;
      for(auto i : nd) vol *= i;
      volCB = parity ? vol : vol / 2;
    }
    MDLatticeParam(const MDLatticeParam<D> &param) :
	   nd{param.Extents()},
     bc(param.GetBC()),
     order(param.GetOrder()),
     parity(param.GetParity()),
     vol(param.GetVol()),
     volCB(param.GetVolCB()) {}

    ~MDLatticeParam() {}

    constexpr int GetNdim(){return Ndim;}
    auto Extents() const {return nd;   }
    auto Extent(const int i) const {return nd[i];}

    BCType    GetBC()   {return bc;   }
    OrderType GetOrder(){return order;}

    inline bool GetParity() const {return parity;}
    inline int  GetVol()    const {return vol;}
    inline int  GetVolCB()  const {return volCB;}
    inline int  GetSurface(const int i, const int j)  const {return nd[i]*nd[j];}

    inline void SetParity(bool parity_){parity = parity_;}

};

template<typename T, int M, int D>
class Grid {
  private :
    const MDLatticeParam<D> &param;
    int nSites;

  public:
    using RegType = T;

    Grid(const MDLatticeParam<D> &param_) : param(param_){
      //number of regs:
      nSites = !param.GetParity() ? param.GetVol() : param.GetVolCB();
    }

    int  NSites() const {return nSites;}
    int  Extent(const int i) const {return param.Extent(i);}
    auto Extents() const {return param.Extents();}
    int  LattExtent(const int i) const {return i == 0 ? param.Extent(i) / M : param.Extent(i);}
    auto LattExtents() const {
      std::array<int, D> latt_dims =  Extents();
      latt_dims[0] /= M;
      return latt_dims;
    }

    ~Grid(){}
};

template<typename T, int M, int D, typename Arg>
class MDLattice : public Grid<T, M, D> {

  using MDLattParam = MDLatticeParam<D>;

  private :

  std::vector<std::array<T, M>> v;
  //Aux fields:
  std::vector<std::array<T, M>> tmp1;
  std::vector<std::array<T, M>> tmp2;
  //Model specific args
  Arg &arg;//arg includes MDLatticeParam.

  public:
    MDLattice(const MDLattParam &param, Arg &arg) : Grid<T, M, D>(param), v(param.GetVol() / M), tmp1(param.GetVol() / M), tmp2(param.GetVol() / M), arg(arg) {}
    MDLattice(const Grid<T, M, D> &grid, Arg &arg) : Grid<T, M, D>(grid), v(grid.NSites() / M), tmp1(grid.NSites() / M), tmp2(grid.NSites() / M), arg(arg) {}

    ~MDLattice() {}

    Arg &GetArgs() const {return arg;}

    template<typename Policy>
    void SetColdLattice(const Policy &p, const T &&val){
      std::array<T, M> t;
      std::fill(t.begin(), t.end(), val);

      std::generate(p, v.begin(), v.end(), [&t]() {return t;});
      std::generate(p, tmp1.begin(), tmp1.end(), [&t]() {return t;});
      std::generate(p, tmp2.begin(), tmp2.end(), [&t]() {return t;});
    }

    //Collection of default accessors:
    inline std::array<T, M>& GetLattPoint(const int i){ return v[i];}
    inline std::array<T, M>& GetTmp1Point(const int i){ return tmp1[i];}
    inline std::array<T, M>& GetTmp2Point(const int i){ return tmp2[i];}

    std::vector<std::array<T, M> >& V()   { return v;   }
    std::vector<std::array<T, M> >& Tmp1(){ return tmp1;}
    std::vector<std::array<T, M> >& Tmp2(){ return tmp2;}
};
