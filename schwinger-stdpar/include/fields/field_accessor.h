#pragma once
//
#include <core/color_spinor.h>

template<GenericFieldTp field_tp, bool is_constant = false, int bSize = 1>
class FieldAccessor {
  public:
    using data_tp  = typename field_tp::data_tp;
    
    static constexpr std::size_t nDir    = field_tp::nDir;
    static constexpr std::size_t nDim    = field_tp::nDim;    
    static constexpr std::size_t nColor  = field_tp::nColor;
    static constexpr std::size_t nSpin   = field_tp::nSpin;    
    static constexpr std::size_t nParity = field_tp::nParity; 
    
    static constexpr FieldType type = field_tp::type();       

    using SpinorTp = impl::ColorSpinor<data_tp, nColor, nSpin, bSize>;
    using LinkTp   = impl::ColorMatrix<data_tp, nColor, bSize>;
    
    using Indices  = std::make_index_sequence<nDim>;      

    using AccessorTp      = typename std::remove_cvref_t< decltype( std::declval<field_tp>().template Accessor<is_constant>()) >;  
    using GhostAccessorTp = typename std::remove_cvref_t< decltype( std::declval<field_tp>().template Accessor<is_constant>()) >;      
  
    AccessorTp       field_accessor;
    GhostAccessorTp  ghost_accessor;    
    
    FieldAccessor(const field_tp &field ) : field_accessor(field.template Accessor<is_constant>()), 
                                            ghost_accessor(field.template Accessor<is_constant>()) {}        
      
    /**
       @brief 2-d accessor functor
       @param[in] x coords
       @param[in] y coords
       @param[in] i spin dof or dir      
       @return Complex number at this spin and color index
    */
    template <int nparity = nParity>
    requires (nparity == 1)
    inline data_tp& operator()(int x, int y, int i) { return field_accessor(x,y,i); }

    /**
       @brief 2-d accessor functor
       @param[in] x coords
       @param[in] y coords
       @param[in] i spin dof or dir       
       @return Complex number at this spin and color index
    */    
    template <int nparity = nParity>
    requires (nparity == 1)    
    inline const data_tp& operator()(int x, int y, int i) const { return field_accessor(x,y,i); }
    
    /**
       @brief 2-d accessor functor
       @param[in] x coords
       @param[in] y coords
       @param[in] i spin dof or dir      
       @return Complex number at this spin and color index
    */
    template <int nparity = nParity>
    requires (nparity == 2)
    inline data_tp& operator()(int x, int y, int i, int p) { return field_accessor(x,y,d,p); }

    /**
       @brief 2-d accessor functor
       @param[in] x coords
       @param[in] y coords
       @param[in] i spin dof or dir      
       @return Complex number at this spin and color index
    */
    template <int nparity = nParity>
    requires (nparity == 2)    
    inline const data_tp& operator()(int x, int y, int i, int p) const { return field_accessor(x,y,d,p); } 
    //
    
    template<std::size_t... Idxs, int nspin = nSpin>
    requires (nspin == 2) 
    inline decltype(auto) load_spinor(std::index_sequence<Idxs...>, const std::array<int, nDim>& x) const {
      std::array<data_tp, nColor*nSpin*bSize> spinor{this->field_accessor(x[Idxs]..., 0), this->field_accessor(x[Idxs]..., 1)};    
      return SpinorTp(spinor);//2-component spinor
    }    
        
    template<std::size_t... Idxs, int ncolor = nColor>
    requires (ncolor == 1) 
    inline decltype(auto) load_parity_link(std::index_sequence<Idxs...>, const std::array<int, nDim>& x, const int &d, const int &parity) const {
      std::array<data_tp, nColor*nColor*bSize> link{this->field_accessor(x[Idxs]..., d, parity)};
      return LinkTp(link);
    }      
    
    template<FieldType field_type = type>
    requires (field_type == FieldType::VectorFieldType)    
    inline decltype(auto) operator()(const std::array<int, nDim> &x, const int &d, const int &p) const {
      return load_parity_link(Indices{}, x, d, p);    
    }            

    template<FieldType field_type = type>
    requires (field_type == FieldType::SpinorFieldType)        
    inline decltype(auto) operator()(const std::array<int, nDim> &x) const {
      return load_spinor(Indices{}, x);    
    }  
    
    inline decltype(auto) Extent(int d) const { return field_accessor.extent(d); }              
};



