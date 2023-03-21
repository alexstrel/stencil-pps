#pragma once

  template <template<GenericContainerTp, typename> class field>
  class BlockField : public std::vector<field> {
    using block_field  = std::vector<field>;
    //
    using container_tp = typename field::container_tp;
    using data_tp      = typename field::data_tp;    

    template <typename U> block_field make_set(std::vector<U> &v) { return block_field{v.begin(), v.end()}; }
    
    template <typename U> block_field make_set(BlockField<U>  &v) { return block_field{v.begin(), v.end()}; }    

  public:
    //
    BlockField()                   = default;
    BlockField(const BlockField &) = default;
    BlockField(BlockField &&)      = default;

    template <typename U> BlockField(U &v) {
      auto vset = make_set(v);
      
      block_field::reserve(vset.size());
      block_field::insert(block_field::end(), vset.begin(), vset.end());
    }

    template <IteratorTp I> BlockField(I first, I last) {
      block_field::reserve(last - first);
      for (auto it = first; it != last; it++) block_field::push_back(*it);
    }

    T& operator[](size_t idx) const { return block_field::operator[](idx).get(); }

  };

  template <template<GenericContainerTp, typename> class field> using CBlockField = const BlockField<field>;


