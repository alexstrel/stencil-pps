#pragma once

using vector_tp         = std::vector<std::complex<Float>>;
using sloppy_vector_tp  = std::vector<std::complex<float>>;

using pmr_vector_tp         = std::pmr::vector<std::complex<Float>>;
using sloppy_pmr_vector_tp  = std::pmr::vector<std::complex<float>>;

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
void run_simple_pmr_dslash(auto &&transformer, auto params, const int X, const int T, const int niter) {
  //
  const auto sf_args = SpinorFieldArgs<nSpin>{{X, T}, {0, 0, 0, 0}};
  //
  const auto gf_args    = GaugeFieldArgs<nDir>{{X, T}, {0, 0}};
  //
  auto gauge            = create_field<vector_tp, decltype(gf_args)>(gf_args);
  //
  init_u1(gauge);
  //using low precision dslash:  
  constexpr bool copy_gauge = true;
  //  
  auto sloppy_gauge = create_field<decltype(gauge), sloppy_vector_tp, copy_gauge>(gauge);  
  //
  auto src_spinor = create_field_with_buffer<pmr_vector_tp, decltype(sf_args)>(sf_args);
  auto dst_spinor = create_field_with_buffer<pmr_vector_tp, decltype(sf_args)>(sf_args);
  auto chk_spinor = create_field_with_buffer<pmr_vector_tp, decltype(sf_args)>(sf_args);
  //
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

  for(int i = 0; i < niter; i++) {
    // Apply dslash	  
    mat(dst_spinor, src_spinor, transformer, FieldOrder::LexFieldOrder);
  }
  
  //DslashRef(chk_spinor, src_spinor, gauge, params.M, params.r, {X, T});  
  
  gauge.destroy();
  sloppy_gauge.destroy();
  
  dst_spinor.destroy();
  src_spinor.destroy();

}

template<int nDir, int nSpin, int N>
void run_mrhs_pmr_dslash(auto &&transformer, auto params, const int X, const int T, const int niter) {
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
  using pmr_spinor_t  = Field<pmr_vector_tp, decltype(sf_args)>;//
  //
  constexpr bool use_pmr_buffer = true;
  //
  auto src_block_spinor = create_block_spinor< pmr_spinor_t, decltype(sf_args), use_pmr_buffer >(sf_args, N); 
  //
  for (int i = 0; i < src_block_spinor.Size(); i++) init_spinor( src_block_spinor.v[i] );
  
  auto dst_block_spinor = create_block_spinor< pmr_spinor_t, decltype(sf_args), use_pmr_buffer >(sf_args, N); 
 
  for(int i = 0; i < niter; i++) {
    mat(dst_block_spinor, src_block_spinor);
  } 

  src_block_spinor.destroy();
  dst_block_spinor.destroy();  
  gauge.destroy();
}

