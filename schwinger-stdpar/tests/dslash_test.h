#pragma once

template<int nDir, int nSpin>
void run_simple_dslash(auto &&transformer, auto params, const int X, const int T, const int niter) {
  //
  const auto sf_args = SpinorFieldArgs<nSpin>{{X, T}, {0, 0}};
  //
  auto src_spinor = create_field<vector_tp, decltype(sf_args)>(sf_args);
  auto dst_spinor = create_field<vector_tp, decltype(sf_args)>(sf_args);
  //auto chk_spinor = create_field<vector_tp, decltype(sf_args)>(sf_args);
  //
  const auto gf_args    = GaugeFieldArgs<nDir>{{X, T}, {0, 0}};
  //
  auto gauge      = create_field<vector_tp, decltype(gf_args)>(gf_args);  
  //
  init_u1(gauge);
  init_spinor(src_spinor);
  // Setup dslash arguments:
  auto &&u_ref        = gauge.View();
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
  
}

template<int nDir, int nSpin>
void run_eo_dslash(auto &&transformer, auto params, const int X, const int T, const int niter) {
  //
  const auto eo_sf_args = SpinorFieldArgs<nSpin>{{X/2, T}, {0, 0}, FieldOrder::EOFieldOrder};
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

  for(int i = 0; i < niter; i++) {
    // Apply dslash	  
    eo_mat(eo_dst_spinor, eo_src_spinor, transformer, FieldOrder::EOFieldOrder);
  }  
}

template<int nDir, int nSpin, int N>
void run_mrhs_dslash(auto &&transformer, auto params, const int X, const int T, const int niter) {
  //
  const auto sf_args = SpinorFieldArgs<nSpin>{{X, T}, {0, 0}};
  //
  //auto src_spinor = create_field<vector_tp, decltype(sf_args)>(sf_args);
  //
  const auto gf_args    = GaugeFieldArgs<nDir>{{X, T}, {0, 0}};
  //
  auto gauge = create_field<vector_tp, decltype(gf_args)>(gf_args);    
  //
  init_u1(gauge);
  //init_spinor(src_spinor);
  // Setup dslash arguments:
  auto &&u_ref        = gauge.View();
  using gauge_tp      = decltype(gauge.View());

  std::unique_ptr<DslashArgs<gauge_tp, nSpin>> dslash_args_ptr(new DslashArgs{u_ref});

  auto &dslash_args = *dslash_args_ptr;

  // Create dslash matrix
  auto mat = Mat<decltype(dslash_args), Dslash, decltype(params)>{dslash_args, params};    
  //
  using spinor_t  = Field<vector_tp, decltype(sf_args)>;//
  //
  auto src_block_spinor = BlockSpinor< spinor_t, decltype(sf_args) >{sf_args, N};
  //
  for (int i = 0; i < src_block_spinor.Size(); i++) init_spinor( src_block_spinor.v[i] );
  
  auto dst_block_spinor = BlockSpinor< spinor_t, decltype(sf_args) >{sf_args, N};

  for(int i = 0; i < niter; i++) {
    mat(dst_block_spinor, src_block_spinor);
  } 
}

