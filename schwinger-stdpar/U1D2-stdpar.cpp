#include <common.h>
#include <field.h>
#include <dslash_factory.h>
//
using Float   = float;

//constexpr int D = 2;
//constexpr int N = 2;

constexpr int ldim = 16384;
constexpr int tdim = 8192;

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()s
std::uniform_real_distribution<Float> dis(0.f, 1.f);

void init_u1(auto &field){
   for (auto &i : field.Data()) i = std::polar(1.0f,dis(gen));	
}

void init_spinor(auto &field){
   for (auto &i : field.Data()) i = std::complex<Float>(1.f, 0.f);
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
  const Float mass = 0.05;
  DslashParam<Float> dslash_param{1.0 / (2.0*(mass +2.0))};

  // allocate and initialize the working lattices, matrices, and vectors
  //
  const auto sf_args = SpinorFieldArgs{ldim, tdim};
  const auto gf_args = GaugeFieldArgs{ldim, tdim};
  //
  auto src_spinor = Field<std::vector<std::complex<Float>>, decltype(sf_args)>(sf_args);
  auto dst_spinor = Field<std::vector<std::complex<Float>>, decltype(sf_args)>(sf_args);
  auto gauge      = Field<std::vector<std::complex<Float>>, decltype(gf_args)>(gf_args);

  init_u1(gauge);
  init_spinor(src_spinor);
  
  // Setup dslash arguments:
  auto &&u_ref    = gauge.Get();
  using gauge_tp  = typename std::remove_cvref_t<decltype(u_ref)>;

  constexpr std::size_t nspin = 2;

  std::unique_ptr<DslashArgs<gauge_tp, nspin>> dslash_args_ptr(new DslashArgs{u_ref});

  auto &dslash_args = *dslash_args_ptr;

  // Create dslash matrix
  auto mat = Mat<Dslash<decltype(dslash_args)>, decltype(dslash_args), decltype(dslash_param)>{dslash_args, dslash_param};

  const int niter = 1000;

  for(int i = 0; i < niter; i++) {
    // Apply dslash	  
    mat(dst_spinor, src_spinor);
  }

  // initialize the data
  bool verbose = true;
  
  if (verbose > 0) {
    std::cout << "Number of sites = " << ldim << " x " << tdim << "." << std::endl;
    std::cout << std::flush;
  }

  return 0;
}
