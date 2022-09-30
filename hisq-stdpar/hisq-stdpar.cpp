#include <common.h>
#include <ranges>
#include <vector>

#include <field.h>

#define LDIM 16
#define TDIM 16

#include <random>

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()s
std::uniform_real_distribution<> dis(-1.f, 1.f);
//
using Float = float;

constexpr int N     = 3;
constexpr int ndirs = 4;
//
constexpr int bSize = 32;
constexpr int Nrows = N;
constexpr int Ncols = N;


//--------------------------------------------------------------------------------
// initializes su3_matrix
template<ArithmeticTp T, int n , int b>
void init_mat(hisq::sun_matrix<T, n, b> &m) {
  T re = dis(gen);
  T im = dis(gen);
#pragma unroll
  for(int k = 0; k < hisq::sun_matrix<T,n,b>::Rows(); ++k){
#pragma unroll  
    for(int l = 0; l < hisq::sun_matrix<T,n,b>::Cols(); ++l) {
#pragma unroll	  
      for(int j = 0; j < b; j++){       
         m(k, l, j) = std::complex<T>(re, im);
      }
    }
  }

  return;
}

//--------------------------------------------------------------------------------
// initializes su3_vector
template<ArithmeticTp T, int n , int b>
void init_vec(hisq::sun_vector<T, n, b> &v) {
  T re = dis(gen);
  T im = dis(gen);
#pragma unroll
  for(int k = 0; k < v.size(); ++k) {
#pragma unroll  
    for(int j = 0; j < b; j++){    
      v(k, j) = std::complex<T>(re, im);
    }
  }
  return;
}

//--------------------------------------------------------------------------------
// initialize lattice data
template<ArithmeticTp T, int n, int b>
void make_data(std::vector<hisq::sun_vector<T, n, b>> &src, std::vector<hisq::sun_matrix<T, n, b>> &fat, std::vector<hisq::sun_matrix<T, n, b>> &lng) {
  for(int i = 0; i < src.size(); i++) {
    init_vec(src[i]);
#pragma unroll    
    for(int dir = 0; dir < ndirs; dir++){
      init_mat(fat[ndirs*i + dir]);
      init_mat(lng[ndirs*i + dir]);
    }
  }
}

//--------------------------------------------------------------------------------
int main(int argc, char **argv)
{

  using su3_vector = hisq::sun_vector<Float, N, bSize>;
  using su3_matrix = hisq::sun_matrix<Float, N, bSize>;  
  
  if(argc < 2){
    std::cerr << "Usage <workgroup size>" << std::endl;
    exit(1);
  }

  int workgroup_size = atoi(argv[1]);//bsize

  int iterations = 100;
  int ldim       = LDIM;
  int tdim       = TDIM;  

  // allocate and initialize the working lattices, matrices, and vectors
  int volume = (ldim*ldim*ldim*tdim) / bSize;
  //
  Field<su3_vector> src(volume);
  Field<su3_vector> dst(volume);
  Field<su3_matrix> fat(volume*ndirs);
  Field<su3_matrix> lng(volume*ndirs);

  // initialize the data
  make_data<Float, N, bSize>(src.Even(), fat.Even(), lng.Even());
  make_data<Float, N, bSize>(src.Odd(), fat.Odd(), lng.Odd());

  bool verbose = true;
  
  if (verbose > 0) {
    std::cout << "Number of sites = " << ldim << "^3 x " << tdim << "." << std::endl;
    std::cout << "Executing " << iterations << " iterations" << std::endl;
    std::cout << std::flush;
  }
#if 0  
  if (verbose > 0)
    std::cout << "Total execution time = " << ttotal << " secs" << std::endl;
  // calculate flops/s, etc.
  // each matrix vector multiply is 3*(12 mult + 12 add) = (36 mult + 36 add) = 72 ops
  // sixteen mat vec operations per site 16*72 = 1152
  // plus 15 complex sums = 1152+30 = 1182 ops per site
  const double tflop = (double)iterations * even_sites_on_node * 1182;
  std::cout << "Total GFLOP/s = " << tflop / ttotal / 1.0e9 << std::endl;

  // calculate ideal, per site, data movement for the dslash kernel
  const double memory_usage = (double)even_sites_on_node *
    (sizeof(su3_matrix)*4*4   // 4 su3_matrix reads, per direction
     +sizeof(su3_vector)*16    // 16 su3_vector reads
     +sizeof(size_t)*16        // 16 indirect address reads
     +sizeof(su3_vector));     // 1 su3_vector write
  std::cout << "Total GByte/s (GPU memory) = " << iterations * memory_usage / ttotal / 1.0e9 << std::endl;

  // check memory usage
  if (verbose > 0) {
    const double memory_allocated = (double)total_sites *
      (sizeof(su3_matrix)*4*4   // 4 su3_matrices, each 4x total_sites in size
       +sizeof(su3_vector)*2     // 2 su3_vectors
       +sizeof(size_t)*4*4);     // 4 index arrays with 4 directional dimensions each
    std::cout << "Total allocation for matrices = " << memory_allocated / 1048576.0 << std::endl;
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0)
      std::cout << "Approximate memory usage = " << usage.ru_maxrss/1024.0 << std::endl;
  }
#endif
  return 0;
}
