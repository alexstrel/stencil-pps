#include <iostream>
#include <cassert>
#include <cmath>
#include <stdlib.h>
#include <iomanip>

#include <time.h>
#include <sys/time.h>

#include <algorithm>
#include <vector>
#include <memory>
#include <numeric>
#include <execution>
#include <chrono>

#include <lattice.h>
#include <stencil_params.h>
#include <stencil_algorithms.h>
#include <stencil_utils.h>

//#define DEBUG_STENCIL
int main(){
  constexpr auto dims  = 3;
  constexpr auto gridsz= 256;
#ifdef __NVCOMPILER_CUDA__
  constexpr auto inner_range = 1;//better 2 for double and 4 for floats
#else
  constexpr auto inner_range = gridsz / 4;
#endif
  constexpr auto stencil_type = StencilType::FaceCentered;

  using Float=float;
  using HKParams = HK3DParams<Float, stencil_type>;

  const std::array<int, dims> nd{gridsz, gridsz, gridsz};

  const MDLatticeParam<dims> param(nd);
  Grid<Float, inner_range, dims> grid(param);

  const Float kappa     = 0.1;
  const Float length    = 1000.0;
  const Float tinterval = 0.5;
  const int   nsteps    = 6553;
  HKParams HK3DArgs(param, kappa, length, tinterval, nsteps);
  //
  MDLattice<Float, inner_range, dims, decltype(HK3DArgs)> hk3d_ref(grid, HK3DArgs);
  //set volumes:
  auto policy = std::execution::par_unseq;
  //hk3d_ref.SetColdLattice(policy, 1.0);

  std::cout << "Start initialization :: " << std::endl;

  create_field<Float, inner_range>(hk3d_ref.Tmp1(), hk3d_ref.Extents(), HK3DArgs.dl, kappa, length, 0.0);

  std::cout << "Done initialization :: " << param.GetVol() << std::endl;

  constexpr StencilPolicy stencil_policy = StencilPolicy::DefaultPolicy;
  const int outer_range = param.GetVol() / inner_range;

  printf("Running forward Euler iterations for the 3d heat kernel %d times with params [c0=%le , c1=%le] ..\n", nsteps, HK3DArgs.c0, HK3DArgs.c1);
  fflush(stdout);

  // Declare timers
  std::chrono::high_resolution_clock::time_point time_begin, time_end;

  FwdEulerIters<Float, inner_range, decltype(policy), stencil_type, stencil_policy, dims, HK3DParams> hk_fwd_euler_algo(policy, HK3DArgs, hk3d_ref.V(), hk3d_ref.Tmp1(), outer_range);

  time_begin = std::chrono::high_resolution_clock::now();
  //launch the compute task
  hk_fwd_euler_algo.apply(nsteps);
  //
  time_end = std::chrono::high_resolution_clock::now();

  Float time = nsteps * HK3DArgs.dt;

  auto &f_final = nsteps & 1 ? hk3d_ref.V() : hk3d_ref.Tmp1();

  create_field<Float, inner_range>(hk3d_ref.Tmp2(), hk3d_ref.Extents(), HK3DArgs.dl, kappa, length, time);

  double err = accuracy<Float, inner_range>(hk3d_ref.Tmp2(),f_final);

  double elapsed_time = std::chrono::duration_cast<std::chrono::duration<double> >(time_end - time_begin).count();

  double Gflops = param.GetVol()*(stencil_type == StencilType::FaceCentered ? 8.0 : stencil_type == StencilType::FaceEdgeCentered ? 21.0 : 30.0)*nsteps/elapsed_time * 1.0e-09;
  double Gstens = param.GetVol()*1.0*nsteps/elapsed_time * 1.0e-06;
  double thput  = param.GetVol() * sizeof(Float) * 3.0 * nsteps
      / elapsed_time * 1.0e-09;

  fprintf(stderr, "Elapsed time : %.3f (s)\n", elapsed_time);
  fprintf(stderr, "FLOPS        : %.3f (GFlops)\n", Gflops);
  fprintf(stderr, "Updates      : %.3f (Mupdates/sec)\n", Gstens);
  fprintf(stderr, "Throughput   : %.3f (GB/s)\n", thput);
  fprintf(stderr, "Accuracy     : %e\n", err);

  return 0;
}
