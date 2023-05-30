#include <common.h>
#include <field.h>
#include <dslash_factory.h>
//
#include <dslash_ref.h>
//
using Float   = float;

//constexpr int D = 2;
//constexpr int N = 2;

constexpr int ldim = 4096;
constexpr int tdim = 4096;

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

   Float tol = 1e-12;
   
   for (int i = 0; i < src1.GetLength(); i++){
//   for (int i = 0; i < 1024; i++){   
     if (norm(src1.Data()[i] - src2.Data()[i]) > tol) 
     {
       std::cout <<  src1.Data()[i] << " :: " << src2.Data()[i] << " " << i << std::endl;
     } 
   }
}

void fill(auto &field_accessor) {
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
  check_spinor(chk_spinor, dst_spinor); 
    
  // initialize the data
  bool verbose = true;
  
  if (verbose > 0) {
    std::cout << "Number of sites = " << ldim << " x " << tdim << "." << std::endl;
    std::cout << std::flush;
  }

  return 0;
}
