#pragma once

enum class StencilType {
  FaceCentered           = 7,  //3d - 7pt  sentcils
  FaceEdgeCentered       = 19, //3d - 19pt stencils
  FaceEdgeCornerCentered = 27  //3d - 27pt stencils
};

enum class Shift  {NoShift=0, ShiftXp1=+1, ShiftXm1=-1, ShiftYp1=+2, ShiftYm1=-2, ShiftZp1=+3, ShiftZm1=-3};

constexpr std::array<Shift, 6> face_shifts{Shift::ShiftXp1, Shift::ShiftXm1, Shift::ShiftYp1, Shift::ShiftYm1, Shift::ShiftZp1, Shift::ShiftZm1};

