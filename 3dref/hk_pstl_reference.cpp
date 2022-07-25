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

#include <ranges>

#include <stencil_utils.h>
#include <stencil_impl.h>

using Float = float;

constexpr auto stencil_type = StencilType::FaceEdgeCornerCentered;

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

int main(){
  Float c0, c1;

  c1 = kappa*dt/(dx*dx);

  if(stencil_type == StencilType::FaceCentered){
    c0 = 1.0 - 6.0*c1;
  } else if (stencil_type == StencilType::FaceEdgeCentered) {
    c0 = 1.0 - 4.0*c1;
    c1 = c1 / 3.0;//?
  } else if (stencil_type == StencilType::FaceEdgeCornerCentered) {
    c0 = 1.0 - (44.0 / 13)*c1;
    c1 = (3.0*c1) / 13.0;//?
  }

  //set execution policy:
  auto policy = std::execution::par_unseq;

  //initialize fields
  std::vector<Float> v1(vol);
  std::vector<Float> v2(vol);
  //
  std::vector<Float> v3(vol);
  //
  create_field<Float>(v3, nd, kappa, length, 0.0);
  //
  std::fill(policy, v2.begin(), v2.end(), 0.0);

  std::copy(policy, v3.begin(), v3.end(), v1.begin());

  std::cout << "Done initialization :: " << vol << std::endl;

  printf("Running forward Euler iterations for the 3d heat kernel %d times with params [c0=%le , c1=%le] ..\n", nsteps, c0, c1);
  fflush(stdout);

  struct timeval time_begin, time_end;

  //Create stencil functor instances:
  std::unique_ptr<Generic3DStencil<Float, stencil_type>> func_ptr(new Generic3DStencil<Float, stencil_type>(c0, c1, nd));

  auto outer_loop_range = std::ranges::views::iota(0, vol);
  
  gettimeofday(&time_begin, NULL);
  const int check_interval = 250;
  //launch iterations
  for(int k = 0; k < nsteps; k++) {
    func_ptr->SetData(v1); 
    if(k % check_interval != 0) {  
      std::for_each(policy,
                    std::ranges::begin(outer_loop_range),
                    std::ranges::end(outer_loop_range),
                    [&stencil= *func_ptr, out = v2.data()] (const int i) {out[i] = stencil(i);});
    } else {
      double l2nrm = std::transform_reduce( policy, 
		      			    begin(v1), 
					    end(v1), 
					    begin(v2), 
					    0., 
					    std::plus<double>(), 
					    [&func= *func_ptr] (const auto &in_el, 
						                      auto &out_el) 
					    {return func(out_el, in_el);} );

      std::cout << "L2 norm (k= " << k <<" ) :: " << sqrt(l2nrm) << std::endl;
    }
    
    std::swap(v1, v2);
  }

  gettimeofday(&time_end, NULL);

  Float time = nsteps * dt;

  auto &f_final = nsteps & 1 ? v2 : v1;
  auto &f_tmp   = nsteps & 1 ? v1 : v2;
  create_field<Float>(f_tmp, nd, kappa, length, time);

  double err = accuracy<Float>(f_tmp,f_final);

  double elapsed_time = (time_end.tv_sec - time_begin.tv_sec)+(time_end.tv_usec - time_begin.tv_usec)*1.0e-6;
  double Gflops = vol*(stencil_type == StencilType::FaceCentered ? 8.0 : stencil_type == StencilType::FaceEdgeCentered ? 21.0 : 30.0)*nsteps/elapsed_time * 1.0e-09;
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
