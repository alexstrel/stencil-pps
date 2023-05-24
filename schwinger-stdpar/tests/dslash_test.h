#pragma once

#include <chrono>

using vector_tp     = std::vector<std::complex<Float>>;

//using sloppy_pmr_vector_tp  = std::pmr::vector<std::complex<float>>;
//using pmr_vector_tp         = std::pmr::vector<std::complex<Float>>;

template<typename Float>
void DslashRef(auto &out_spinor, const auto &in_spinor, const auto &gauge_field, const Float mass, const Float r, const std::array<int, 2> n) {//const int nx, const int ny
  
  const Float constant = (mass + 2.0*r);
  
  const int nx = n[0];
  const int ny = n[1];
  //  
  std::complex<Float> tmp = std::complex<Float>(0.);
  
  auto I = [](auto x){ return Float(-x.imag(), x.real());};  
  
  auto out          = out_spinor.Accessor();
  const auto in     = in_spinor.Accessor();
  const auto gauge  = gauge_field.Accessor();  
  
  for(int x = 0; x < nx; x++) {
    const int xp1 = (x+1) == nx ? 0    : (x+1);
    const int xm1 = (x-1) == -1 ? nx-1 : (x-1);
    
    for(int y = 0; y < ny; y++) {
      const int yp1 = (y+1) == ny ? 0    : (y+1);
      const int ym1 = (y-1) == -1 ? ny-1 : (y-1);

      const Float fwd_bndr = yp1 == 0 ? -1.0 : 1.0;
      const Float bwd_bndr = y   == 0 ? -1.0 : 1.0;  

      //
      tmp = constant * in(x,y,0)

	- 0.5*(gauge(x,y,0) * (in(xp1,y,0) - in(xp1,y,1)) + conj(gauge(xm1,y,0)) * (in(xm1,y,0) + in(xm1,y,1)))
	
	- 0.5*(gauge(x,y,1) * fwd_bndr*(in(x,yp1,0) + I(in(x,yp1,1))) + conj(gauge(x,ym1,1)) * bwd_bndr*(in(x,ym1,0) - I(in(x,ym1,1))));
      
      out(x,y,0) = tmp;
      
      //
      tmp = constant * in(x,y,1) 

	- 0.5*(gauge(x,y,0) * (in(xp1,y,1) - in(xp1,y,0)) + conj(gauge(xm1,y,0)) * (in(xm1,y,1) + in(xm1,y,0)))
	
	- 0.5*(gauge(x,y,1) * fwd_bndr*(in(x,yp1,1) - I(in(x,yp1,0))) + conj(gauge(x,ym1,1)) * bwd_bndr*(in(x,ym1,1) + I(in(x,ym1,0))));
      
      out(x,y,1) = tmp;
    }
  }
}

template<int nDir, int nSpin>
void run_dslash_test(auto params, const int X, const int T, const int niter) {
  //
  constexpr int nSpinorParity = 1;
  constexpr int nGaugeParity  = 2;
  //  
  const auto cs_param = SpinorFieldArgs<nSpinorParity>{{X/2, T}, {0, 0, 0, 0}, FieldParity::EvenFieldParity};
  //
  auto src_spinor = create_field<vector_tp, decltype(cs_param)>(cs_param);
  auto dst_spinor = create_field<vector_tp, decltype(cs_param)>(cs_param);
  //
  const auto gauge_param = GaugeFieldArgs<nGaugeParity>{{X, T}, {0, 0}}; 
  //
  auto gauge = create_field<vector_tp, decltype(gauge_param)>(gauge_param);
  //
  init_u1(gauge);
  init_spinor(src_spinor);  

  auto &&u_ref   = gauge.View();
  using gauge_tp = decltype(gauge.View());
    
  std::unique_ptr<DslashArgs<gauge_tp, nSpin>> dslash_args_ptr(new DslashArgs{u_ref});

  auto &dslash_args = *dslash_args_ptr;

  // Create dslash matrix
  auto mat = Mat<decltype(dslash_args), Dslash>{dslash_args};
  //
  const auto s1 = params.M + static_cast<Float>(2.0)*params.r;
  const auto s2 = static_cast<Float>(0.5);

  auto transformer = [=](const auto &x, const auto &y) {return (s1*x-s2*y);};      

  auto wall_start = std::chrono::high_resolution_clock::now();

  for(int i = 0; i < niter; i++) {
    // Apply dslash	  
    mat(dst_spinor, src_spinor, transformer, FieldParity::EvenFieldParity);
  } 
 
  auto wall_stop = std::chrono::high_resolution_clock::now();

  auto wall_diff = wall_stop - wall_start;
  
  auto wall_time = (static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(wall_diff).count()) / 1e6)  / niter;

  std::cout << "Done for EO version : time per iteration is > " << wall_time << "sec." << std::endl;

  gauge.destroy();
  
  dst_spinor.destroy();
  src_spinor.destroy();   
}

template<int nDir, int nSpin, int N, bool use_pmr_buffer = false>
void run_mrhs_dslash_test(auto params, const int X, const int T, const int niter) {
  //
  constexpr int nSpinorParity = 1;
  constexpr int nGaugeParity  = 2;  
  //  
  const auto cs_param    = SpinorFieldArgs<nSpinorParity>{{X/2, T}, {0, 0, 0, 0}, FieldParity::EvenFieldParity};
  //
  const auto gauge_param = GaugeFieldArgs<nGaugeParity>{{X, T}, {0, 0}};
  //
  auto gauge = create_field<vector_tp, decltype(gauge_param)>(gauge_param);    
  //  
  init_u1(gauge);
  //
  auto &&u_ref        = gauge.View();
  using gauge_tp      = decltype(gauge.View());

  std::unique_ptr<DslashArgs<gauge_tp, nSpin>> dslash_args_ptr(new DslashArgs{u_ref});

  auto &dslash_args = *dslash_args_ptr;

  // Create dslash matrix
  auto mat = Mat<decltype(dslash_args), Dslash>{dslash_args};    
  //
  using spinor_t  = Field<vector_tp, decltype(cs_param)>;//
  //
  auto src_block_spinor = create_block_spinor< spinor_t, decltype(cs_param), use_pmr_buffer>(cs_param, N); 
  //
  for (int i = 0; i < src_block_spinor.nComponents(); i++) init_spinor( src_block_spinor.v[i] );
  
  auto dst_block_spinor = create_block_spinor< spinor_t, decltype(cs_param), use_pmr_buffer>(cs_param, N); 

  const auto s1 = params.M + static_cast<Float>(2.0)*params.r;
  const auto s2 = static_cast<Float>(0.5);

  auto transformer = [=](const auto &x, const auto &y) {return (s1*x-s2*y);};

  auto wall_start = std::chrono::high_resolution_clock::now(); 

  for(int i = 0; i < niter; i++) {
    mat(dst_block_spinor, src_block_spinor, transformer, FieldParity::EvenFieldParity);
  } 

  auto wall_stop = std::chrono::high_resolution_clock::now();

  auto wall_diff = wall_stop - wall_start;

  auto wall_time = (static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(wall_diff).count()) / 1e6)  / niter;

  std::cout << "Done for MRHS version (N =  " << N << ") : time per iteration is > " << wall_time << "sec." << std::endl;

  src_block_spinor.destroy();
  dst_block_spinor.destroy();  
  gauge.destroy();
}

