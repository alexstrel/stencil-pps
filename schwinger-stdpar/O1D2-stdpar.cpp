#include <common.h>
#include <field.h>
//
using Float   = float;

//constexpr int D = 2;
//constexpr int N = 2;

constexpr int ldim = 16;
constexpr int tdim = 16;

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()s
std::uniform_real_distribution<Float> dis(0.f, 1.f);

void init_u1(auto &field){
   for (auto &i : field.Get()) i = std::polar(1.0f,dis(gen));	
}

void init_spinor(auto &field){
   for (auto &i : field.Get()) i = std::complex<Float>(1.f, 0.f);
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
   std::cout << "Print components for field : " << field.Get().data() << std::endl;

   auto print = [](const auto& e) { std::cout << "Element " << e << std::endl; };

   std::for_each(field.Get().begin(), field.Get().begin()+range, print);
}

//--------------------------------------------------------------------------------
int main(int argc, char **argv)
{
  // allocate and initialize the working lattices, matrices, and vectors
  //
  const SpinorFieldArgs sf_args{ldim, tdim};
  const auto gf_args = GaugeFieldArgs{ldim, tdim};
  //
  auto src_spinor = Field<std::vector<std::complex<Float>>, decltype(sf_args)>(sf_args);
  auto dst_spinor = Field<std::vector<std::complex<Float>>, decltype(sf_args)>(sf_args);
  auto gauge      = Field<std::vector<std::complex<Float>>, decltype(gf_args)>(gf_args);

  init_u1(gauge);
  init_spinor(src_spinor);

  print_range(gauge, 4);

  auto odd_gauge = gauge.Odd();
  print_range(odd_gauge, 4);

  auto even_gauge= gauge.Even();
  print_range(even_gauge, 4); 

  auto even_gauge_acc = gauge.Even().Accessor();

  std::cout << even_gauge_acc(2,0,0) << std::endl;
  //
  fill(even_gauge_acc);
  //
  std::cout << even_gauge_acc(2,0,0) << std::endl;

  // initialize the data
  bool verbose = true;
  
  if (verbose > 0) {
    std::cout << "Number of sites = " << ldim << " x " << tdim << "." << std::endl;
    std::cout << std::flush;
  }

  return 0;
}
