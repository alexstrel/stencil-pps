#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>
#include <oneapi/dpl/iterator>
#include <oneapi/dpl/random>

#include <iostream>
#include <cassert>
#include <cmath>
#include <stdlib.h>
#include <iomanip>

#include <time.h>
#include <sys/time.h>

#include <vector>
#include <memory>
#include <numeric>
#include <limits>
#include <chrono>

#include <CL/sycl.hpp>
using oneapi::dpl::counting_iterator;

#include <stencil_utils.h>
#include <stencil_impl.h>

using Float = float;

constexpr auto stencil_type = StencilTp::FaceEdgeCornerCentered;

const int dims  = 3;
const int gridsz= 128;

const std::array<int, dims> nd{gridsz, gridsz, gridsz};
const int vol = nd[0]*nd[1]*nd[2];

const Float kappa     = 0.1f;
const Float length    = 1000.0f;
const Float tinterval = 0.5f;
const int   nsteps    = 100;
const Float dt        = tinterval / nsteps;
const Float dx        = length / (nd[0]+static_cast<Float>(1.0f));

template <bool enable_gpu_backend = false>
decltype(auto) get_queue(){

  if constexpr ( enable_gpu_backend == false) {
    std::cout << "WARNING: clang++ generated x86 backend. For Intel GPUs use -DUSE_GPU option.\n" << std::endl;
    return sycl::queue(sycl::cpu_selector{});
  } else {
    std::cout << "WARNING: clang++ generated GPU backend.\n" << std::endl;
    return sycl::queue(sycl::gpu_selector{});
  }
}

template <typename T, typename queue_t>
decltype(auto) get_allocator(queue_t& cq){
  return sycl::usm_allocator<T, sycl::usm::alloc::shared>{cq};
}

template<typename lambda_tp, typename policy_tp>
void dispatch_stencil_kernel(const lambda_tp& site_stencil_kernel, const int vol, const policy_tp &policy){
  //	
  auto outer_loop_range = vol;
  //
  std::for_each(policy,
                counting_iterator(0),
                counting_iterator(outer_loop_range),
                site_stencil_kernel);

  return;
}

int main(){
  
  Float c0, c1, c2, c3;

  c1 = kappa*dt/(dx*dx);

  if(stencil_type == StencilTp::FaceCentered){
    c0 = 1.0 - 6.0*c1;
  } else if (stencil_type == StencilTp::FaceEdgeCentered) {
    c0 = 1.0 - 4.0*c1;
    c1 = c1 / 3.0;//?
  } else if (stencil_type == StencilTp::FaceEdgeCornerCentered) {
    c0 = 1.0 - (44.0 / 13.0)*c1;
    c1 = (3.0*c1) / 13.0;//?
  }
  // extra:
  c2 = c1 * 0.5;
  c3 = c1 / 3.0;

  //set execution policy:
  auto cq = get_queue<true>();
  //
  auto policy = oneapi::dpl::execution::make_device_policy(cq);
  //
  auto alloc = get_allocator<Float, decltype(cq)>(cq);
  //
  std::vector<Float, decltype(alloc)> v1(vol, alloc);
  std::vector<Float, decltype(alloc)> v2(vol, alloc);
  //
  //initialize fields:
  create_field<decltype(policy), Float, decltype(alloc)>(policy, v1, nd, kappa, length, 0.f);
  //
  std::fill(policy, v2.begin(), v2.end(), 0.f);

  std::cout << "Done initialization :: " << vol << std::endl;

  printf("Running forward Euler iterations for the 3d heat kernel %d times with params [c0=%le , c1=%le, c2=%le, c3=%le] ..\n", nsteps, c0, c1, c2, c3);
  fflush(stdout);

  struct timeval time_begin, time_end;
  
  using StencilArgs = GenericNDStencilArg<Float, decltype(alloc), dims, Float, Float, Float, Float>;
  
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
  auto stencil_kernel = [stencil = func_ptr.get()] (const int i) { stencil->operator()(i); };

  for(int k = 0; k < nsteps; k++) {
    // 
    if((k+1) % check_interval != 0) {  
      dispatch_stencil_kernel(stencil_kernel, vol, policy);
    } else {    //??
 #if 0
      double l2nrm = std::transform_reduce( policy, 
                                            begin(v1), 
                                            end(v1), 
                                            begin(v2), 
                                            0.,  
                                            std::plus<double>{}, 
                                            [&func= *func_ptr] (const auto &in_el, auto &out_el) 
					    			       {return func(out_el, in_el);} );

      std::cout << "L2 norm (k= " << k <<" ) :: " << sqrt(l2nrm) << std::endl;
#endif
    }
    
    args.Swap();
  }

  gettimeofday(&time_end, NULL);

  Float time = nsteps * dt;

  auto &f_final = nsteps & 1 ? v2 : v1;
  auto &f_tmp   = nsteps & 1 ? v1 : v2;
  create_field<decltype(policy), Float, decltype(alloc)>(policy, f_tmp, nd, kappa, length, time);

  double err = accuracy<Float, decltype(alloc)>(f_tmp,f_final);

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
