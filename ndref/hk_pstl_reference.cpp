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
#include <limits>
#include <execution>
#include <chrono>

#include <ranges>

#include <stencil_utils.h>
#include <stencil_impl.h>

constexpr double _1div3_ = 1.0 / 3.0;

using Float = float;

constexpr auto stencil_type = StencilTp::FaceEdgeCornerCentered;

const int dims  = 3;
const int gridsz= 768;

const std::array<int, dims> nd{gridsz, gridsz, gridsz};
const int vol = nd[0]*nd[1]*nd[2];

const Float kappa     = 0.1;
const Float length    = 1000.0;
const Float tinterval = 0.5;
const int   nsteps    = 100;
const Float dt        = tinterval / nsteps;
const Float dx        = length / (nd[0]+static_cast<Float>(1.0));

void dispatch_stencil_kernel(auto&& site_stencil_kernel, const int vol){
  //	
  auto policy = std::execution::par_unseq;
  //
  auto outer_loop_range = std::ranges::views::iota(0, vol);
  //
  std::for_each(policy,
                std::ranges::begin(outer_loop_range),
                std::ranges::end(outer_loop_range),
                site_stencil_kernel);

  return;
}

int main(){

  using StencilArgs = GenericNDStencilArg<Float, dims, Float, Float, Float, Float>;//NB!
  
  Float c0, c1, c2, c3;

  c1 = kappa*dt/(dx*dx);

  if(stencil_type == StencilTp::FaceCentered){
    c0 = 1.0 - 6.0*c1;
  } else if (stencil_type == StencilTp::FaceEdgeCentered) {
    c0 = 1.0 - 4.0*c1;
    c1 = c1 / 3.0;//?
  } else if (stencil_type == StencilTp::FaceEdgeCornerCentered) {
    c0 = 1.0 - (44.0 / 13)*c1;
    c1 = (3.0*c1) / 13.0;//?
  }
  // extra:
  c2 = c1 * 0.5;
  c3 = c1 *_1div3_;
  
  //set execution policy:
  auto policy = std::execution::par_unseq;

  std::vector<Float> v1(vol);
  std::vector<Float> v2(vol);
  //
    //initialize fields:
  create_field<decltype(policy), Float>(policy, v1, nd, kappa, length, 0.0);
  //
  std::fill(policy, v2.begin(), v2.end(), 0.0);

  std::cout << "Done initialization :: " << vol << std::endl;

  printf("Running forward Euler iterations for the 3d heat kernel %d times with params [c0=%le , c1=%le] ..\n", nsteps, c0, c1);
  fflush(stdout);

  struct timeval time_begin, time_end;
  
  //Create stencil functor instances:
  std::unique_ptr<StencilArgs> args_ptr(new StencilArgs{v2, v1, nd, c0, c1, c2, c3});
  
  auto& args = *args_ptr;
  
  auto func_ptr = std::make_shared<GenericNDStencil<stencil_type, StencilArgs>>(args);
  
  gettimeofday(&time_begin, NULL);
#ifdef __NVCOMPILER_CUDA__
  const int check_interval = 250;
#else
  const int check_interval = std::numeric_limits<int>::max();
#endif  
  //launch iterations
  auto stencil_kernel = [&stencil = *func_ptr] (const auto i) { stencil(i); };

  for(int k = 0; k < nsteps; k++) {
    // 
    if((k+1) % check_interval != 0) {  
      dispatch_stencil_kernel(stencil_kernel, vol);
    } else {    
      double l2nrm = std::transform_reduce( policy, 
                                            begin(v1), 
                                            end(v1), 
                                            begin(v2), 
                                            0.,  
                                            std::plus<double>{}, 
                                            [&func= *func_ptr] (const auto &in_el, auto &out_el) 
					    			       {return func(out_el, in_el);} );

      std::cout << "L2 norm (k= " << k <<" ) :: " << sqrt(l2nrm) << std::endl;
    }
    
    args.Swap();
  }

  gettimeofday(&time_end, NULL);

  Float time = nsteps * dt;

  auto &f_final = nsteps & 1 ? v2 : v1;
  auto &f_tmp   = nsteps & 1 ? v1 : v2;
  create_field<decltype(policy), Float>(policy, f_tmp, nd, kappa, length, time);

  double err = accuracy<Float>(f_tmp,f_final);

  double elapsed_time = (time_end.tv_sec - time_begin.tv_sec)+(time_end.tv_usec - time_begin.tv_usec)*1.0e-6;
  double Gflops = vol*(stencil_type == StencilTp::FaceCentered ? 8.0 : stencil_type == StencilTp::FaceEdgeCentered ? 21.0 : 30.0)*nsteps/elapsed_time * 1.0e-09;
  double Gstens = vol*1.0*nsteps/elapsed_time * 1.0e-06;
  double thput  = vol * sizeof(Float) * 3.0 * nsteps
      / elapsed_time * 1.0e-09;

  fprintf(stderr, "Elapsed time : %.3f (s)\n", elapsed_time);
  fprintf(stderr, "FLOPS        : %.3f (GFlops)\n", Gflops);
  fprintf(stderr, "Updates      : %.3f (Mupdates/sec)\n", Gstens);
  fprintf(stderr, "Throughput   : %.3f (GB/s)\n", thput);
  fprintf(stderr, "Accuracy     : %e\n", err);

  return 0;
}
