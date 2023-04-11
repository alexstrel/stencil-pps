#pragma once
#include <common.h>
#include <memory>
#include <memory_resource>

static std::size_t pool_bytes = 0ul;

static std::shared_ptr<std::byte[]> pool_ptr = nullptr;

/**
*      @return mapped memory allocated
*         
**/
//std::size_t allocated();

/**
 *     @return is prefetching support enabled
 **/
//bool is_prefetch_enabled();


inline void allocate_pmr_pool(const std::size_t bytes) {
  //
  pool_bytes = bytes;
  //
  if(pool_ptr == nullptr) pool_ptr = std::make_shared<std::byte[]>(pool_bytes);  
} 

inline void release_pmr_pool() {
  pool_ptr.reset();

  pool_ptr  = nullptr;
  pool_bytes= 0ul;   
}

decltype(auto) get_pmr_pool(const size_t bytes, const std::size_t offset = 0ul) {
  if(pool_ptr == nullptr) {
    std::cerr << "Cannot create a pmr resource, upstream pool does not exist." << std::endl;
    exit(-1);
  }

  if (bytes > pool_bytes) {
    std::cerr << "Cannot create a pmr resource, requested number of pool bytes to big." << std::endl;
    exit(-1); 
  }
 
  return std::pmr::monotonic_buffer_resource{pool_ptr.get() + offset, bytes};
}

template<PMRContainerTp T>
decltype(auto) get_pmr_container(const size_t n) {
  if(pool_ptr == nullptr) {
    std::cerr << "Cannot create a pmr resource, upstream pool does not exist." << std::endl;
    exit(-1);
  }

  const std::size_t bytes = n*sizeof(typename T::value_type);

  if (bytes > pool_bytes) {
    std::cerr << "Cannot create a pmr resource, requested number of pool bytes to big." << std::endl;
    exit(-1);
  }

  auto pmr_pool = std::pmr::monotonic_buffer_resource{pool_ptr.get(), bytes};

  return T(n, &pmr_pool);
}



