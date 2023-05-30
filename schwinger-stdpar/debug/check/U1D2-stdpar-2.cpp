#include <common.h>
#include <field.h>
#include <dslash_factory.h>
//
#include <dslash_ref.h>
//
using Float   = float;

//constexpr int D = 2;
//constexpr int N = 2;

constexpr int ldim = 32;
constexpr int tdim = 32;

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()s
std::uniform_real_distribution<Float> dis(0.f, 1.f);

void init_u1(auto &field){
   for (auto &i : field.Data()) i = std::polar(1.0f,dis(gen));	
}

void init_spinor(auto &field){
   for (auto &i : field.Data()) i = std::complex<Float>(dis(gen), dis(gen));
}

void check_spinor(auto &src1, auto &src2){

   Float tol = 1e-13;
   
   for (int i = 0; i < src1.GetLength(); i++){
//   for (int i = 0; i < 1024; i++){   
     if (norm(src1.Data()[i] - src2.Data()[i]) > tol) 
     {
       std::cout <<  src1.Data()[i] << " :: " << src2.Data()[i] << " " << i << std::endl;
     } 
   }
}

void fill(auto &src_field_accessor) {
  const int mu = field_accessor.extent(2);
  for(int i = 0; i < field_accessor.extent(0); i++){
    for(int j = 0; j < field_accessor.extent(1); j++){	 
#pragma unroll 
      for(int k = 0; k < mu; k++){
	 field_accessor(i,j,k) = std::complex<Float>(1.0, 0.0);     
      }	       
    }
  }
}

void convert_lex2eo(auto &dst_field_accessor, const auto &src_field_accessor){

  const int mu = src_field_accessor.extent(2);
  const int V  = src_field_accessor.extent(0) * src_field_accessor.extent(1); 
  const int Vh = dst_field_accessor.extent(0) * dst_field_accessor.extent(1); 
  for(int p = 0; p < dst_field_accessor.extent(3); p++){    
    for(int i = 0; i < dst_field_accessor.extent(0); i++){
      for(int j = 0; j < dst_field_accessor.extent(1); j++){
        int parity_bit = j & 1; 	 
        int parity_offset = parity_bit == 0 ? p : (1-p);
#pragma unroll 
        for(int k = 0; k < mu; k++){
	  dst_field_accessor(i,j,k,p) = src_field_accessor(2*i+parity_offset,j,k);     
        }	       
      }
    }
  }
  return;
}

void convert_eo2lex(auto &dst_field_accessor, const auto &src_field_accessor){

  const int mu = src_field_accessor.extent(2);
  const int Vh = src_field_accessor.extent(0) * src_field_accessor.extent(1); 
  const int V  = dst_field_accessor.extent(0) * dst_field_accessor.extent(1); 
  for(int p = 0; p < src_field_accessor.extent(3); p++){    
    for(int i = 0; i < src_field_accessor.extent(0); i++){
      for(int j = 0; j < src_field_accessor.extent(1); j++){
        int parity_bit = j & 1; 	 
        int parity_offset = parity_bit == 0 ? p : (1-p);
#pragma unroll 
        for(int k = 0; k < mu; k++){
	  dst_field_accessor(2*i+parity_offset,j,k) = src_field_accessor(i,j,k,p);     
        }	       
      }
    }
  }
  return;
}

void print_range(auto &field, const int range){
   std::cout << "Print components for field : " << field.Data().data() << std::endl;

   auto print = [](const auto& e) { std::cout << "Element " << e << std::endl; };

   std::for_each(field.Data().begin(), field.Data().begin()+range, print);
}

//--------------------------------------------------------------------------------
int main(int argc, char **argv)
{
  //
  using vector_tp = std::vector<std::complex<Float>>;

  const Float mass = 0.05;
  const Float r    = 1.0;

  DslashParam<Float> dslash_param{mass, r};

  constexpr int nDir  = 2;  
  constexpr int nSpin = 2;

  // allocate and initialize the working lattices, matrices, and vectors
  //
  const auto sf_args = SpinorFieldArgs<nSpin>{ldim, tdim};
  //
  auto src_spinor = create_field<vector_tp, decltype(sf_args)>(sf_args);
  auto dst_spinor = create_field<vector_tp, decltype(sf_args)>(sf_args);
  auto chk_spinor = create_field<vector_tp, decltype(sf_args)>(sf_args);  
  
  const auto gf_args = GaugeFieldArgs<nDir>{ldim, tdim};
  //
  auto gauge      = create_field<vector_tp, decltype(gf_args)>(gf_args);

  init_u1(gauge);
  init_spinor(src_spinor);
  
  // Setup dslash arguments:
  auto &&u_ref    = gauge.Reference();
  using gauge_tp  = typename std::remove_cvref_t<decltype(u_ref)>;

  constexpr std::size_t nspin = 2;

  std::unique_ptr<DslashArgs<gauge_tp, nspin>> dslash_args_ptr(new DslashArgs{u_ref});

  auto &dslash_args = *dslash_args_ptr;

  // Create dslash matrix
  auto mat = Mat<decltype(dslash_args), Dslash, decltype(dslash_param)>{dslash_args, dslash_param};  

  const int niter = 1000;
  
  const auto scale1 = mass + static_cast<Float>(2.0)*r;
  const auto scale2 = static_cast<Float>(0.5);

  auto transformer = [=](const auto &x, const auto &y) {return (scale1*x-scale2*y);};  

  for(int i = 0; i < niter; i++) {
    // Apply dslash	  
    mat(dst_spinor, src_spinor, transformer);
  }

  //
  auto &&in_ref    = src_spinor.Reference();
  auto &&out_ref   = dst_spinor.Reference();
  auto &&chk_ref   = chk_spinor.Reference();  
  
  auto out__ = out_ref.Accessor(); 
  auto chk__ = chk_ref.Accessor();     
  
  std::cout << " >>> " << out__(0,0,0) << std::endl;
  //
  Dpsi(chk_ref, in_ref, u_ref, mass, ldim, tdim);
  
  std::cout << " <<< " << chk__(0,0,0) << std::endl;
  //check_spinor(chk_spinor, dst_spinor); 
  
  //now lets do parity spinors:
  const auto eo_gf_args = GaugeFieldArgs<nDir>{ldim, tdim, FieldOrder::EOFieldOrder}; 
  const auto eo_sf_args = SpinorFieldArgs<nSpin>{ldim, tdim, FieldOrder::EOFieldOrder};
  //
  auto eo_gauge  = create_field<vector_tp, decltype(eo_gf_args)>(eo_gf_args);  
  //
  auto eo_src_spinor = create_field<vector_tp, decltype(eo_sf_args)>(eo_sf_args);    
  auto eo_dst_spinor = create_field<vector_tp, decltype(eo_sf_args)>(eo_sf_args);
  //
  //init_spinor(eo_src_spinor); 
  auto&& eo_src_spinor_ref = eo_src_spinor.Reference();  
  auto lll = eo_src_spinor_ref.ExtAccessor();
  auto nnn = in_ref.Accessor();
  convert_lex2eo(lll, nnn);     
  //
  auto [e_src_spinor_ref, o_src_spinor_ref] = eo_src_spinor.EODecompose();    
  //
  auto &&eo_u_ref = eo_gauge.Reference(); 
  //
  auto eo_u  = eo_u_ref.ExtAccessor();   
  auto lex_u = u_ref.Accessor();     
  //
  convert_lex2eo(eo_u, lex_u);   
  //
  std::cout << "<<< ??? <<< " << eo_u(1,1,0,0) << std::endl;  
  std::cout << "<<< ??? <<< " << lex_u(3,1,0) << std::endl;    
  //
  auto ll1 = e_src_spinor_ref.Accessor();
  auto ll2 = in_ref.Accessor();  
  std::cout << "<<< !!ll1 <<< " << ll1(1,1,0) << std::endl;  
  std::cout << "<<< !!ll2 <<< " << ll2(3,1,0) << std::endl;      
  //  
  auto&& e_dst_spinor_ref = eo_dst_spinor.Even();
  //
  Dpsi_parity(e_dst_spinor_ref, o_src_spinor_ref, e_src_spinor_ref, eo_u_ref, mass, ldim / 2, tdim, 0);  
  //
  auto chk2__ = e_dst_spinor_ref.Accessor(); 
  std::cout << " <<< " << chk2__(0,0,0) << std::endl;  
  // initialize the data
  bool verbose = true;
  
  if (verbose > 0) {
    std::cout << "Number of sites = " << ldim << " x " << tdim << "." << std::endl;
    std::cout << std::flush;
  }

  return 0;
}
