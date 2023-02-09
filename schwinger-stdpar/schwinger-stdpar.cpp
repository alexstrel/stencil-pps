#include <common.h>
#include <ranges>
#include <vector>

#include <field.h>

constexpr int D = 2;
constexpr int N = 2;

constexpr int ldim = 16;
constexpr int tdim = 16;

#include <random>

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()s
std::uniform_real_distribution<> dis(0.f, 1.f);
//
using Float   = float;

void init_u1(auto &field){
   for (auto &i : field.Get()) i = std::polar(1.0,dis);	
}

void init_spinor(auto &field){
   for (auto &i : field.Get()) i = std::complex<Float>(1.0, 0.0);
}

//--------------------------------------------------------------------------------
int main(int argc, char **argv)
{
  // allocate and initialize the working lattices, matrices, and vectors
  //
  const auto sf_args = FieldArgs<N, 1>{ldim, tdim};
  const auto gf_args = FieldArgs<1, D>{ldim, tdim};
  //
  Field<std::complex<Float>> src_spinor(sf_args);
  Field<std::complex<Float>> dst_spinor(sf_args);
  Field<std::complex<Float>> gauge(gf_args);

  init_u1(gauge);
  init_spinor(src_spinor);

  // initialize the data
  bool verbose = true;
  
  if (verbose > 0) {
    std::cout << "Number of sites = " << ldim << " x " << tdim << "." << std::endl;
    std::cout << std::flush;
  }

  return 0;
}
