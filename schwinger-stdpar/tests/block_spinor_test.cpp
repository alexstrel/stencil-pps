#include <fields/field.h>
#include <fields/block_field.h>
#include <kernels/dslash_factory.h>

//#include <dslash_test.h>
//
using Float   = double;
//
//using vector_tp = std::vector<std::complex<Float>>;

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()s
std::uniform_real_distribution<Float> dis(0.f, 1.f);

void init_u1(auto &field){
   for (auto &i : field.Data()) i = std::polar(1.0f,dis(gen));	
}

void init_spinor(auto &field){
   for (auto &i : field.Data()) i = dis(gen); //std::complex<Float>(1.f, 0.f);
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

//template<FieldTp field_tp>
template<GenericSpinorFieldTp field_tp>
void print_range(field_tp &field, const int range){
   std::cout << "Print components for field : " << field.Data().data() << " ::: " << field.GetLength() << std::endl;

   auto print = [](const auto& e) { std::cout << "Element " << e <<  "\t address " << &e << std::endl; };
 
   if(range > field.GetLength()) {
     std::cout << "WARNING: print range exceeds field length, nop." << std::endl;
     return;
   }

   std::for_each(field.Data().begin(), field.Data().begin()+range, print);
}

template<int N, bool to_eo = true>
void convert_field(auto &dst_field, auto &src_field){
   const auto [Nx, Ny] = src_field.GetDims();

   auto X = std::views::iota(0, Nx);
   auto Y = std::views::iota(0, Ny);

   auto idx = std::views::cartesian_product(Y, X);//Y is the slowest index, X is the fastest

   auto &&dst = dst_field.View();
   auto &&src = src_field.View(); 

   std::for_each(std::execution::par_unseq,
                 idx.begin(),
                 idx.end(), [=](const auto i) {
		      auto [y, x] = i;
		      const int parity = (x + y) & 1;
		      // 
		      auto dstU = dst.ExtAccessor();
		      auto srcU = src.template Accessor<to_eo>();
#pragma unroll
		      for(int j = 0; j < N; j++){
		        if constexpr (to_eo) {  
		          dstU(x/2,y,j,parity) = srcU(x,y,j);
			} else {
			  srcU(x,y,j) = dstU(x/2,y,j,parity);
	                }		
		      }
		   });

}

#include <memory.h>

//--------------------------------------------------------------------------------
int main(int argc, char **argv)
{
  //
  constexpr int nDir  = 2;  
  constexpr int nSpin = 2;  
  constexpr int N     = 2;

  constexpr int X = 1024;
  constexpr int T = 1024;

  const int niter = 1000;

  using vector_tp = std::vector<std::complex<Float>>;

  const auto sargs = SpinorFieldArgs<nSpin>{{X, T}, {0, 0, 0, 0}};
  //
  auto src_spinor = create_field<vector_tp, decltype(sargs)>(sargs);  

  using pmr_vector_tp = std::pmr::vector<std::complex<Float>>;
  
  std::cout << "Create PMR spinor" << std::endl;
  //
  auto pmr_src_spinor = create_field_with_buffer<pmr_vector_tp, decltype(sargs)>(sargs);
  //
  pmr_src_spinor.show();
  //
  init_spinor(pmr_src_spinor);
  //
  print_range(pmr_src_spinor, 4);
  //
  using pmr_vector_lp_tp = std::pmr::vector<std::complex<float>>;

  std::cout << "Exported PMR spinor..." << std::endl;

  constexpr bool destroy_src = false;

  auto next_pmr_spinor = export_pmr_field<decltype(pmr_src_spinor), pmr_vector_lp_tp>(pmr_src_spinor, destroy_src);  
  next_pmr_spinor.show();

  print_range(next_pmr_spinor, 4);

  std::cout << "NOW they share same buffer!" << std::endl;

  //pmr_src_spinor.destroy();
  //
  pmr_src_spinor.show();
  //
  print_range(pmr_src_spinor, 4);

  std::cout << "Create Regular Block spinor" << std::endl;
  //
  using spinor_t = Field<vector_tp,  decltype(sargs)>;
  //
  auto block_src_spinor = create_block_spinor< spinor_t, decltype(sargs)>(sargs, N);
  //
  block_src_spinor[0].show();
  //
  for(int i = 0; i < block_src_spinor.Size(); i++) {
    init_spinor(block_src_spinor[i]);
    print_range(block_src_spinor[i], 4);
    //
    block_src_spinor[i].show();
  }

  std::cout << "Create PMR Block spinor" << std::endl;
  //
  using pmr_spinor_t = Field<pmr_vector_tp,  decltype(sargs)>;

  auto pmr_block_src_spinor = create_block_spinor< pmr_spinor_t, decltype(sargs), true>(sargs, N);
  //
  pmr_block_src_spinor[0].show();

  for(int i = 0; i < pmr_block_src_spinor.Size(); i++) {
    init_spinor(pmr_block_src_spinor[i]);
    print_range(pmr_block_src_spinor[i], 4);
    //
    pmr_block_src_spinor[i].show();
  }

  print_range(pmr_block_src_spinor[0], 4);

  using arg_tp = typename std::remove_cvref_t<decltype(pmr_block_src_spinor.ExportArg())>;

  auto next_pmr_block_src_spinor = export_pmr_block_spinor< decltype(pmr_block_src_spinor) >(pmr_block_src_spinor, N);

//  auto new_pmr_block_spinor = create_block_spinor_with_buffer< pmr_spinor_t, decltype(pmr_sargs) >(pmr_sargs, N);
{
  std::cout << "NEW block field data:: " << std::endl;
  next_pmr_block_src_spinor[0].show();

  for(int i = 0; i < next_pmr_block_src_spinor.Size(); i++) {
    print_range(next_pmr_block_src_spinor[i], 4);
    //
    next_pmr_block_src_spinor[i].show();
  }
}
  //
  pmr_block_src_spinor[0].show();

  std::cout << "Now check block field data:: " << std::endl;
  for(int i = 0; i < pmr_block_src_spinor.Size(); i++) {
    print_range(pmr_block_src_spinor[i], 4);
    //
    pmr_block_src_spinor[i].show();
  }
    //

  // initialize the data
  bool verbose = true;
  
  if (verbose > 0) {
    std::cout << "Number of sites = " << X << " x " << T << "." << std::endl;
    std::cout << std::flush;
  }

  return 0;
}
