#pragma once

template <ReferenceFieldTp ref_field>
class block_field_ref : public std::vector<field> {
    using block_field  = std::vector<ref_field>;
    //
    using container_tp = typename ref_field::container_tp;
    using data_tp      = typename ref_field::data_tp;    

    template <AllocatedFieldTp alc_field_tp> block_field make_set(std::vector<alc_field_tp> &v) { 
      using ref_field_tp = decltype(std::declval<alc_field_tp>().Reference());
          
      return block_field{v.begin(), v.end()}; }
    
    template <ReferenceFieldTp ref_field_tp> block_field make_set(BlockField<ref_field_tp>  &v) { return block_field{v.begin(), v.end()}; }    

  public:
    //
    BlockField()                   = default;
    BlockField(const BlockField &) = default;
    BlockField(BlockField &&)      = default;

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


