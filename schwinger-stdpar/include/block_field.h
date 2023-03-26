#pragma once

template <GenericSpinorFieldTp spinor_field>
class block_spinor_field : public std::vector<spinor_field> {
    using block_field  = std::vector<spinor_field>;
    //
    using container_tp = typename spinor_field::container_tp;
    using data_tp      = typename spinor_field::data_tp;    

    template <SpinorFieldTp field_tp> block_field make_set(std::vector<field_tp> &v)         { return block_field{v.begin(), v.end()}; }
    
    template <SpinorFieldTp field_tp> block_field make_set(block_spinor_field<field_tp>  &v) { return block_field{v.begin(), v.end()}; }    

  public:
    //
    block_spinor_field()                           = default;
    block_spinor_field(const block_spinor_field &) = default;
    block_spinor_field(block_spinor_field &&)      = default;

    template <AllocatedFieldTp alloc_field_tp> 
    BlockField(std::vector<alloc_field_tp> &v) {
      using ref_field_tp = decltype(std::declval<alloc_field_tp>().Reference());
      
      block_field<ref_field_tp> vset;
      //
      vset.reserve(v.size());
      
      auto vset = make_set(v);
      vector::reserve(vset.size());
            
      for (auto &f : v) vset.push_back(f.Reference());
    }

    template <IteratorTp I> BlockField(I first, I last) {
      block_field::reserve(last - first);
      for (auto it = first; it != last; it++) block_field::push_back(*it);
    }

    T& operator[](size_t idx) const { return block_field::operator[](idx).get(); }
};

  template <typename field> using CBlockField = const BlockField<field>;


