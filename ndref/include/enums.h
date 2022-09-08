#pragma once

enum class StencilTp {
  FaceCentered           = 7,  //3d - 7pt  sentcils
  FaceEdgeCentered       = 19, //3d - 19pt stencils
  FaceEdgeCornerCentered = 27  //3d - 27pt stencils
};

enum class Shift  {NoShift=0, ShiftXp1=+1, ShiftXm1=-1, ShiftYp1=+2, ShiftYm1=-2, ShiftZp1=+3, ShiftZm1=-3};

enum class Face   {NoFace=0, FaceXp=1, FaceXm=2, FaceYp=4, FaceYm=8, FaceZp=16, FaceZm=32};

enum class KernelType {Interior_Kernel=0, Exterior_Kernel_X=1, Exterior_Kernel_Y=2, Exterior_Kernel_Z=4, Exterior_Kernel_All=7};

constexpr std::array<Shift, 6> shifts{Shift::ShiftXp1, Shift::ShiftXm1, Shift::ShiftYp1, Shift::ShiftYm1, Shift::ShiftZp1, Shift::ShiftZm1};

constexpr std::array<Face, 6>  faces{Face::FaceXp, Face::FaceXm, Face::FaceYp, Face::FaceYm, Face::FaceZp, Face::FaceZm};


