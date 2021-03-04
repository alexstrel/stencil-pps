#pragma once

template<typename T, StencilType StTp>
struct HK3DParams{
  const MDLatticeParam<3> &lattice;

  T c0;
  T c1;

  const std::array<T, 3> dl;//dl[0]=dx,dl[1]=dy,dl[2]=dz
  const T dt;

  HK3DParams(const MDLatticeParam<3> &l, const T kappa, const T length, const T tinterval, const int nsteps) : lattice(l), dl{length / (l.Extent(0)+1.0), length / (l.Extent(1)+1.0), length / (l.Extent(2)+1.0)}, dt( tinterval / nsteps) {
    c1 = kappa*dt/(dl[0]*dl[0]);
    if(StTp == StencilType::FaceCentered){
      c0 = 1.0 - 6*c1;
    } else if (StTp == StencilType::FaceEdgeCentered) {
      c0 = 1.0 - 4*c1;
      c1 = c1 / 3.0;//?
    } else if (StTp == StencilType::FaceEdgeCornerCentered) {
      c0 = 1.0 - (44.0 / 13)*c1;
      c1 = (3.0*c1) / 13.0;//?
    }
    return;
  }
};
