#pragma once

#include <enums.h>
#include <iterators.h>
#include <stencil_params.h>
#include <stencil_impl.h>

//requires c++20 flag
//template <typename RegTp, typename Policy, StencilType StTp, StencilPolicy StPl, int ND, template<typename data_t, StencilType stencil_t> class Args, std::enable_if<std::is_floating_point<RegTp>::value, bool> = true>
template <typename RegTp, int M, typename Policy, StencilType StTp, StencilPolicy StPl, int ND, template<typename data_t, StencilType stencil_t> class Args>
class FwdEulerIters{
  private:
    Policy &policy;

    const Args<RegTp, StTp> &args;

    std::vector<std::array<RegTp, M>> &v1;
    std::vector<std::array<RegTp, M>> &v2;

    const int outer_range;

  public:
    FwdEulerIters(Policy &policy, const Args<RegTp, StTp> &args, std::vector<std::array<RegTp, M>> &f1, std::vector<std::array<RegTp, M>> &f2, const int outer_range) :
      policy(policy),
      args(args),
      v1(f1),
      v2(f2),
      outer_range(outer_range){  }

    void apply(const int nsteps){
      //Create stencil functor instances:
      std::unique_ptr<GenericNDStencil<RegTp, M, StTp, ND>> even_t_func_ptr(new GenericNDStencil<RegTp, M, StTp, ND>(v1, v2, args.c0, args.c1, args.lattice.Extents()));
      std::unique_ptr<GenericNDStencil<RegTp, M, StTp, ND>> odd_t_func_ptr (new GenericNDStencil<RegTp, M, StTp, ND>(v2, v1, args.c0, args.c1, args.lattice.Extents()));
      //launch iterations
      for(int k = 0; k < nsteps; k++) {
        auto &func = (k & 1) == 0 ? *even_t_func_ptr : *odd_t_func_ptr;

        std::for_each(policy,
                      impl::counting_iterator(0),
                      impl::counting_iterator(outer_range),
                      [&func] (const int j) {return func.template operator()<StPl>(j);});
      }
      return;
    }

    ~FwdEulerIters(){}
};
