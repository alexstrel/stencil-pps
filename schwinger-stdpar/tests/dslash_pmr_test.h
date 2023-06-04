#pragma once

using vector_tp         = std::vector<std::complex<Float>>;
using sloppy_vector_tp  = std::vector<std::complex<float>>;

using pmr_vector_tp         = impl::pmr::vector<std::complex<Float>>;
using sloppy_pmr_vector_tp  = impl::pmr::vector<std::complex<float>>;

template<typename Float>
void DslashRef(auto &out_spinor, const auto &in_spinor, const auto &accum_spinor, const auto &gauge_field, const Float mass, const Float r, const std::array<int, 2> n, const int parity) {//const int nx, const int ny
  
  const Float constant = (mass + 2.0*r);
  
  const int nxh = n[0];
  const int ny  = n[1];
  //  
  std::complex<Float> tmp = std::complex<Float>(0.);
  
  auto I = [](auto x){ return Float(-x.imag(), x.real());};  
  
  auto out          = out_spinor.ParityAccessor();
  const auto in     = in_spinor.ParityAccessor();
  //
  const auto gauge  = gauge_field.Accessor(); 
  
  const int other_parity = 1 - parity; 
  
  for(int y = 0; y < ny; y++) {
    const int yp1 = (y+1) == ny ? 0    : (y+1);
    const int ym1 = (y-1) == -1 ? ny-1 : (y-1);

    const Float fwd_bndr = yp1 == 0 ? -1.0 : 1.0;
    const Float bwd_bndr = y   == 0 ? -1.0 : 1.0;
    
    const int parity_bit = y & 1; 
  
    const int fwd_stride = parity_bit ? +1 :  0; 
    const int bwd_stride = parity_bit ?  0 : +1;       
  
    for(int x = 0; x < nxh; x++) {
      //      
      const int xp1 = (x+fwd_stride) == nxh ? 0     : (x+fwd_stride);
      const int xm1 = (x-bwd_stride) == -1  ? nxh-1 : (x-bwd_stride);      
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


void run_pmr_dslash_test(auto params, const int X, const int T, const int niter) {
  //
  constexpr int nSpinorParity = 1;
  constexpr int nGaugeParity  = 2;
  //  
  const auto cs_param = SpinorFieldArgs<nSpinorParity>{{X/2, T}, {0, 0, 0, 0}, FieldParity::EvenFieldParity};
  //
  const auto gauge_param = GaugeFieldArgs<nGaugeParity>{{X, T}, {0, 0}};
  //
  auto gauge             = create_field<vector_tp, decltype(gauge_param)>(gauge_param);
  //
  init_u1(gauge);
  //using low precision dslash:  
  constexpr bool copy_gauge = true;
  //  
  auto sloppy_gauge = create_field<decltype(gauge), sloppy_vector_tp, copy_gauge>(gauge);  
  //
  auto src_spinor  = create_field_with_buffer<sloppy_pmr_vector_tp, decltype(cs_param)>(cs_param);
  auto accum_spinor= create_field_with_buffer<sloppy_pmr_vector_tp, decltype(cs_param)>(cs_param);  
  auto dst_spinor  = create_field_with_buffer<sloppy_pmr_vector_tp, decltype(cs_param)>(cs_param);
  auto chk_spinor  = create_field_with_buffer<sloppy_pmr_vector_tp, decltype(cs_param)>(cs_param);  
  //
  init_spinor(src_spinor);  
  init_spinor(accum_spinor);    
  
  // Setup dslash arguments:
  auto &&u_ref        = sloppy_gauge.View();
  //u_ref.destroy();
  using sloppy_gauge_tp = decltype(sloppy_gauge.View());
  using spinor_ref_tp   = decltype(src_spinor.View());

  std::unique_ptr<DslashArgs<sloppy_gauge_tp, decltype(cs_param)::nspin>> dslash_args_ptr(new DslashArgs{u_ref});

  auto &dslash_args = *dslash_args_ptr;

  // Create dslash matrix
  auto mat = Mat<decltype(dslash_args), Dslash>{dslash_args};   
  
  const auto s1 = static_cast<float>(params.M) + static_cast<float>(2.0)*static_cast<float>(params.r);
  const auto s2 = static_cast<float>(0.5);

  auto transformer = [=](const auto &x, const auto &y) {return (s1*x-s2*y);};
  //
  //mat(dst_spinor, src_spinor, accum_spinor, transformer, FieldParity::EvenFieldParity);        

  auto wall_start = std::chrono::high_resolution_clock::now();   

  for(int i = 0; i < niter; i++) {
    // Apply dslash	  
    mat(dst_spinor, src_spinor, accum_spinor, transformer, FieldParity::EvenFieldParity);
  }
  
  auto wall_stop = std::chrono::high_resolution_clock::now();

  auto wall_diff = wall_stop - wall_start;
  
  auto wall_time = (static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(wall_diff).count()) / 1e6)  / niter;

  std::cout << "Done for EO version : time per iteration is > " << wall_time << "sec." << std::endl; 
  
  src_spinor.show();
  
  dst_spinor.destroy();
  src_spinor.destroy();
  accum_spinor.destroy();  

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto src_spinor2   = create_field_with_buffer<sloppy_pmr_vector_tp, decltype(cs_param)>(cs_param);
  auto dst_spinor2   = create_field_with_buffer<sloppy_pmr_vector_tp, decltype(cs_param)>(cs_param);
  auto accum_spinor2 = create_field_with_buffer<sloppy_pmr_vector_tp, decltype(cs_param)>(cs_param);  
  
  auto wall_start2 = std::chrono::high_resolution_clock::now();   

  for(int i = 0; i < niter; i++) {
    // Apply dslash	  
    mat(dst_spinor2, src_spinor2, accum_spinor2, transformer, FieldParity::EvenFieldParity);
  }
  
  auto wall_stop2 = std::chrono::high_resolution_clock::now();

  auto wall_diff2 = wall_stop2 - wall_start2;
  
  auto wall_time2 = (static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(wall_diff2).count()) / 1e6)  / niter;

  std::cout << "Done for another EO version : time per iteration is > " << wall_time2 << "sec." << std::endl;  

  gauge.destroy();
  sloppy_gauge.destroy();
  
  src_spinor.show();  
  
  dst_spinor2.destroy();
  src_spinor2.destroy();  
  accum_spinor2.destroy();    
}

template<int N>
void run_mrhs_pmr_dslash_test(auto params, const int X, const int T, const int niter) {
  //
  constexpr int nSpinorParity = 1;
  constexpr int nGaugeParity  = 2;
  //  
  const auto cs_param = SpinorFieldArgs<nSpinorParity>{{X/2, T}, {0, 0, 0, 0}, FieldParity::EvenFieldParity};
  //
  const auto gauge_param = GaugeFieldArgs<nGaugeParity>{{X, T}, {0, 0}};
  //
  auto gauge = create_field<vector_tp, decltype(gauge_param)>(gauge_param);    
  //  
  init_u1(gauge);
  //
  auto &&u_ref        = gauge.View();
  using gauge_tp      = decltype(gauge.View());

  std::unique_ptr<DslashArgs<gauge_tp, decltype(cs_param)::nspin>> dslash_args_ptr(new DslashArgs{u_ref});

  auto &dslash_args = *dslash_args_ptr;

  // Create dslash matrix
  auto mat = Mat<decltype(dslash_args), Dslash>{dslash_args};    
  //
  using pmr_spinor_t  = Field<pmr_vector_tp, decltype(cs_param)>;//
  //
  constexpr bool use_pmr_buffer = true;
  //
  auto src_block_spinor = create_block_spinor< pmr_spinor_t, decltype(cs_param), use_pmr_buffer >(cs_param, N); 
  //
  auto accum_block_spinor = create_block_spinor< pmr_spinor_t, decltype(cs_param), use_pmr_buffer >(cs_param, N);   
  //
  for (int i = 0; i < src_block_spinor.nComponents(); i++) init_spinor( src_block_spinor.v[i] );
  for (int i = 0; i < src_block_spinor.nComponents(); i++) init_spinor( accum_block_spinor.v[i] );  
  
  auto dst_block_spinor = create_block_spinor< pmr_spinor_t, decltype(cs_param), use_pmr_buffer >(cs_param, N); 
  
  const auto s1 = params.M + static_cast<Float>(2.0)*params.r;
  const auto s2 = static_cast<Float>(0.5);

  auto transformer = [=](const auto &x, const auto &y) {return (s1*x-s2*y);}; 
  //
  mat(dst_block_spinor, src_block_spinor, accum_block_spinor, transformer, FieldParity::EvenFieldParity);  

  auto wall_start = std::chrono::high_resolution_clock::now();   
 
  for(int i = 0; i < niter; i++) {
    mat(dst_block_spinor, src_block_spinor, accum_block_spinor, transformer, FieldParity::EvenFieldParity);
  } 

  auto wall_stop = std::chrono::high_resolution_clock::now();

  auto wall_diff = wall_stop - wall_start;

  auto wall_time = (static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(wall_diff).count()) / 1e6)  / niter;

  std::cout << "Done for MRHS version (N =  " << N << ") : time per iteration is > " << wall_time << "sec." << std::endl;

  src_block_spinor.destroy();
  dst_block_spinor.destroy();  
  accum_block_spinor.destroy();
  gauge.destroy();  
}

