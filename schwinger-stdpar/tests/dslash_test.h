#pragma once

#include <chrono>

using vector_tp = std::vector<std::complex<Float>>;

//FLOPS per cite : [2*(6 + (6 + 6 + 6 + 6) + 2) + (6 + 6 + 2)]

template<typename Float>
void DslashRef(auto &out_spinor, const auto &in_spinor, const auto &accum_spinor, const auto &gauge_field, const Float mass, const Float r, const std::array<int, 2> n, const int parity) {//const int nx, const int ny
  
  constexpr bool is_constant = true;
  
  const Float constant = (mass + 2.0*r);
  
  const int nxh = n[0];
  const int ny  = n[1];
  //  
  std::complex<Float> tmp = std::complex<Float>(0.);
  std::complex<Float> hlf = std::complex<Float>(0.5, 0.);  
  
  auto I = [](auto x){ return std::complex<Float>(-x.imag(), x.real());};  
  
  auto out          = out_spinor.ParityAccessor();
  const auto in     = in_spinor.template ParityAccessor<is_constant>();
  const auto accum  = accum_spinor.template ParityAccessor<is_constant>();  
  //
  const auto gauge  = gauge_field.template Accessor<is_constant>(); 

  const int other_parity = 1 - parity; 
  
  for(int y = 0; y < ny; y++) {
    const int yp1 = (y+1) == ny ? 0    : (y+1);
    const int ym1 = (y-1) == -1 ? ny-1 : (y-1);

    const std::complex<Float> fwd_bndr = yp1 == 0 ? std::complex<Float>{-1.0, 0.0} : std::complex<Float>{+1.0, 0.0};
    const std::complex<Float> bwd_bndr = y   == 0 ? std::complex<Float>{-1.0, 0.0} : std::complex<Float>{+1.0, 0.0};
    
    const int parity_bit = y & 1; 
  
    const int fwd_stride = parity_bit ? +1 :  0; 
    const int bwd_stride = parity_bit ?  0 : +1;       
  
    for(int x = 0; x < nxh; x++) {
      //      
      const int xp1 = (x+fwd_stride) == nxh ? 0     : (x+fwd_stride);
      const int xm1 = (x-bwd_stride) == -1  ? nxh-1 : (x-bwd_stride);      
      //
      tmp = constant * accum(x,y,0)

	- hlf*(gauge(x,y,0, parity) * (in(xp1,y,0) - in(xp1,y,1)) + conj(gauge(xm1,y,0, other_parity)) * (in(xm1,y,0) + in(xm1,y,1)))
	
	- hlf*(gauge(x,y,1, parity) * fwd_bndr*(in(x,yp1,0) + I(in(x,yp1,1))) + conj(gauge(x,ym1,1, other_parity)) * bwd_bndr*(in(x,ym1,0) - I(in(x,ym1,1))));
      
      out(x,y,0) = tmp;
      
      //
      tmp = constant * accum(x,y,1) 

	- hlf*(gauge(x,y,0, parity) * (in(xp1,y,1) - in(xp1,y,0)) + conj(gauge(xm1,y,0, other_parity)) * (in(xm1,y,1) + in(xm1,y,0)))
	
	- hlf*(gauge(x,y,1, parity) * fwd_bndr*(in(x,yp1,1) - I(in(x,yp1,0))) + conj(gauge(x,ym1,1, other_parity)) * bwd_bndr*(in(x,ym1,1) + I(in(x,ym1,0))));
      
      out(x,y,1) = tmp;
    }
  }
}

void run_dslash_test(auto params, const int X, const int T, const int niter) {
  //
  constexpr int nSpinorParity = 1;
  constexpr int nGaugeParity  = 2;
  //  
  const auto cs_param = SpinorFieldArgs<nSpinorParity>{{X/2, T}, {0, 0, 0, 0}, FieldParity::EvenFieldParity};
  //
  auto src_spinor   = create_field<vector_tp, decltype(cs_param)>(cs_param);
  auto accum_spinor = create_field<vector_tp, decltype(cs_param)>(cs_param);  
  auto dst_spinor = create_field<vector_tp, decltype(cs_param)>(cs_param);
  auto chk_spinor = create_field<vector_tp, decltype(cs_param)>(cs_param);  
  //
  const auto gauge_param = GaugeFieldArgs<nGaugeParity>{{X, T}, {0, 0}}; 
  //
  auto gauge = create_field<vector_tp, decltype(gauge_param)>(gauge_param);
  //
  init_u1(gauge);
  init_spinor(src_spinor);  
  init_spinor(accum_spinor);    

  auto &&u_ref   = gauge.View();
  using gauge_tp = decltype(gauge.View());
    
  std::unique_ptr<DslashArgs<gauge_tp, decltype(cs_param)::nspin>> dslash_args_ptr(new DslashArgs{u_ref});

  auto &dslash_args = *dslash_args_ptr;

  // Create dslash matrix
  auto mat = Mat<decltype(dslash_args), Dslash>{dslash_args};
  //
  const auto s1 = params.M + static_cast<Float>(2.0)*params.r;
  const auto s2 = static_cast<Float>(0.5);

  auto transformer = [=](const auto &x, const auto &y) {return (s1*x-s2*y);};
  //
  mat(dst_spinor, src_spinor, accum_spinor, transformer, FieldParity::EvenFieldParity);        

  auto wall_start = std::chrono::high_resolution_clock::now();

  for(int i = 0; i < niter; i++) {
    // Apply dslash	  
    mat(dst_spinor, src_spinor, accum_spinor, transformer, FieldParity::EvenFieldParity);
  } 
 
  auto wall_stop = std::chrono::high_resolution_clock::now();

  auto wall_diff = wall_stop - wall_start;
  
  auto wall_time = (static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(wall_diff).count()) / 1e6)  / niter;
  
  DslashRef<Float>(chk_spinor, src_spinor, accum_spinor, gauge, params.M, params.r, {X/2, T}, 0); 

  {
    //check_field(chk_spinor, dst_spinor);
    auto &&chk_ = chk_spinor.ParityAccessor();
    auto &&dst_ = dst_spinor.ParityAccessor();     
    //
    check_field(chk_, dst_, 1e-6);
  }
  
  std::cout << "Done for EO version : time per iteration is > " << wall_time << "sec." << std::endl;

  gauge.destroy();
  
  chk_spinor.destroy();  
  dst_spinor.destroy();
  src_spinor.destroy();   
  accum_spinor.destroy();
}

template<int N, bool use_pmr_buffer = false>
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

  std::unique_ptr<DslashArgs<gauge_tp, decltype(cs_param)::nspin>> dslash_args_ptr(new DslashArgs{u_ref});

  auto &dslash_args = *dslash_args_ptr;

  // Create dslash matrix
  auto mat = Mat<decltype(dslash_args), Dslash>{dslash_args};    
  //
  using spinor_t  = Field<vector_tp, decltype(cs_param)>;//
  //
  auto src_block_spinor   = create_block_spinor< spinor_t, decltype(cs_param), use_pmr_buffer>(cs_param, N); 
  auto accum_block_spinor = create_block_spinor< spinor_t, decltype(cs_param), use_pmr_buffer>(cs_param, N);   
  //
  for (int i = 0; i < src_block_spinor.nComponents(); i++) init_spinor( src_block_spinor.v[i] );
  for (int i = 0; i < src_block_spinor.nComponents(); i++) init_spinor( accum_block_spinor.v[i] );  
  
  auto dst_block_spinor = create_block_spinor< spinor_t, decltype(cs_param), use_pmr_buffer>(cs_param, N); 

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

