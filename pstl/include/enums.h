#pragma once

enum class BCType        {Antiperiodic, Dirichlet, Neumann, Periodic, Undefined};
enum class OrderType     {Lexicographic, EvenOdd, Undefined};
enum class RegisterType  {Float, Float2, Float4, Float8, Float16};
enum class Precision     {Quad, Double, Single, Half};
enum class ExecSpace     {X86, Cuda, Sycl, Kokkos};
enum class KernelType    {Internal, ExternalX, ExternalY, ExternalZ};

//enum class StencilType   {Stencil7pt=7, Stencil19pt=19, Stencil27pt=27};
enum class StencilType {
  FaceCentered = 3579, //1d 3pt stencils, 2d 5pt stencils, 3d 7pt sentcils, 4d 9pt stencils
  FaceCornerCentered = 9, //1d - NA, 2d - 9pt stencils, 3d - NA
  FaceEdgeCentered = 19, //1d - NA, 2d - NA, 3d - 19pt stencils
  FaceEdgeCornerCentered = 27//1d - NA, 2d - NA, 3d - 27pt stencils
};
enum class StencilPolicy {DefaultPolicy, Cache2DBlockingPolicy, Cache2p5DBlockingPolicy, Cache3DBlockingPolicy, Cache3p5DBlockingPolicy};
enum class Shift         {NoShift=0, ShiftXp1=+1, ShiftXm1=-1, ShiftYp1=+2, ShiftYm1=-2, ShiftZp1=+3, ShiftZm1=-3, ShiftTp1=+4, ShiftTm1=-4};
