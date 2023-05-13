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
void run_simple_dslash(auto &&transformer, auto params, const int X, const int T, const int niter) {
  //
  const auto sf_args = SpinorFieldArgs<nSpin>{{X, T}, {0, 0, 0, 0}};
  //
  auto src_spinor = create_field<vector_tp, decltype(sf_args)>(sf_args);
  auto dst_spinor = create_field<vector_tp, decltype(sf_args)>(sf_args);
  auto chk_spinor = create_field<vector_tp, decltype(sf_args)>(sf_args);
  //
  const auto gf_args    = GaugeFieldArgs<nDir>{{X, T}, {0, 0}};
  //
  auto gauge            = create_field<vector_tp, decltype(gf_args)>(gf_args);  
  //
  init_u1(gauge);
  init_spinor(src_spinor);
  // Setup dslash arguments:
  auto &&u_ref        = gauge.View();
  //u_ref.destroy();
  using gauge_tp      = decltype(gauge.View());
  using spinor_ref_tp = decltype(src_spinor.View());

  std::unique_ptr<DslashArgs<gauge_tp, nSpin>> dslash_args_ptr(new DslashArgs{u_ref});

  auto &dslash_args = *dslash_args_ptr;

  // Create dslash matrix
  auto mat = Mat<decltype(dslash_args), Dslash, decltype(params)>{dslash_args, params};    
 
  auto wall_start = std::chrono::high_resolution_clock::now();

  for(int i = 0; i < niter; i++) {
    // Apply dslash	  
    mat(dst_spinor, src_spinor, transformer, FieldOrder::LexFieldOrder);
  }

  auto wall_stop = std::chrono::high_resolution_clock::now();
  //
  auto wall_diff = wall_stop - wall_start;
 
  auto wall_time = (static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(wall_diff).count()) / 1e6)  / niter;

  std::cout << "Done for simple version : time per iteration is > " << wall_time << "sec." << std::endl;
 
  //DslashRef(chk_spinor, src_spinor, gauge, params.M, params.r, {X, T});  
  
  gauge.destroy();
  
  dst_spinor.destroy();
  src_spinor.destroy();

}

template<int nDir, int nSpin>
void run_eo_dslash(auto &&transformer, auto params, const int X, const int T, const int niter) {
  //
  const auto eo_sf_args = SpinorFieldArgs<nSpin>{{X/2, T}, {0, 0, 0, 0}, FieldOrder::EOFieldOrder};
  //
  auto eo_src_spinor = create_field<vector_tp, decltype(eo_sf_args)>(eo_sf_args);
  auto eo_dst_spinor = create_field<vector_tp, decltype(eo_sf_args)>(eo_sf_args);
  //auto eo_chk_spinor = create_field<vector_tp, decltype(eo_sf_args)>(eo_sf_args);
  //
  const auto eo_gf_args = GaugeFieldArgs<nDir>{{X/2, T}, {0, 0}, FieldOrder::EOFieldOrder}; 
  //
  auto eo_gauge = create_field<vector_tp, decltype(eo_gf_args)>(eo_gf_args);
  //
  init_u1(eo_gauge);
  init_spinor(eo_src_spinor);  

  auto &&eo_u_ref   = eo_gauge.View();
  using eo_gauge_tp = decltype(eo_gauge.View());
    
  std::unique_ptr<DslashArgs<eo_gauge_tp, nSpin>> eo_dslash_args_ptr(new DslashArgs{eo_u_ref});

  auto &eo_dslash_args = *eo_dslash_args_ptr;

  // Create dslash matrix
  auto eo_mat = Mat<decltype(eo_dslash_args), Dslash, decltype(params)>{eo_dslash_args, params};  

  auto wall_start = std::chrono::high_resolution_clock::now();

  for(int i = 0; i < niter; i++) {
    // Apply dslash	  
    eo_mat(eo_dst_spinor, eo_src_spinor, transformer, FieldOrder::EOFieldOrder);
  } 
 
  auto wall_stop = std::chrono::high_resolution_clock::now();

  auto wall_diff = wall_stop - wall_start;
  
  auto wall_time = (static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(wall_diff).count()) / 1e6)  / niter;

  std::cout << "Done for EO version : time per iteration is > " << wall_time << "sec." << std::endl;

  eo_gauge.destroy();
  
  eo_dst_spinor.destroy();
  eo_src_spinor.destroy();   
}

template<int nDir, int nSpin, int N, bool use_pmr_buffer = false>
void run_mrhs_dslash(auto &&transformer, auto params, const int X, const int T, const int niter) {
  //
  const auto sf_args = SpinorFieldArgs<nSpin>{{X, T}, {0, 0, 0, 0}};
  //
  const auto gf_args = GaugeFieldArgs<nDir>{{X, T}, {0, 0}};
  //
  auto gauge = create_field<vector_tp, decltype(gf_args)>(gf_args);    
  //  
  init_u1(gauge);
  //
  auto &&u_ref        = gauge.View();
  using gauge_tp      = decltype(gauge.View());

  std::unique_ptr<DslashArgs<gauge_tp, nSpin>> dslash_args_ptr(new DslashArgs{u_ref});

  auto &dslash_args = *dslash_args_ptr;

  // Create dslash matrix
  auto mat = Mat<decltype(dslash_args), Dslash, decltype(params)>{dslash_args, params};    
  //
  using spinor_t  = Field<vector_tp, decltype(sf_args)>;//
  //
  auto src_block_spinor = create_block_spinor< spinor_t, decltype(sf_args), use_pmr_buffer>(sf_args, N); 
  //
  for (int i = 0; i < src_block_spinor.Size(); i++) init_spinor( src_block_spinor.v[i] );
  
  auto dst_block_spinor = create_block_spinor< spinor_t, decltype(sf_args), use_pmr_buffer>(sf_args, N); 

  auto wall_start = std::chrono::high_resolution_clock::now(); 

  for(int i = 0; i < niter; i++) {
    mat(dst_block_spinor, src_block_spinor);
  } 

  auto wall_stop = std::chrono::high_resolution_clock::now();

  auto wall_diff = wall_stop - wall_start;

  auto wall_time = (static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(wall_diff).count()) / 1e6)  / niter;

  std::cout << "Done for MRHS version (N =  " << N << ") : time per iteration is > " << wall_time << "sec." << std::endl;

  src_block_spinor.destroy();
  dst_block_spinor.destroy();  
  gauge.destroy();
}

