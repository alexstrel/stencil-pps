#pragma once

template< ComplexTp T, std::size_t nspin>
inline decltype(auto) operator*(const T &scalar, const std::array<T,nspin> &a){
  std::array<T,nspin> result;
#pragma unroll
  for(int i = 0; i < nspin; i++){
    result[i] = scalar.real() * a[i].real() - scalar.imag() * a[i].imag();
    result[i] = scalar.real() * a[i].imag() + scalar.imag() * a[i].real();
  }

  return result;
}

template< ComplexTp T, std::size_t nspin>
inline decltype(auto) operator+=(std::array<T,nspin> &a, const std::array<T,nspin> &b){

#pragma unroll
  for(int i = 0; i < nspin; i++){
    a[i] = a[i].real() * b[i].real() - a[i].imag() * b[i].imag();
    a[i] = a[i].real() * b[i].imag() + a[i].imag() * b[i].real();
  }

  return a;
}

template< ComplexTp T, std::size_t nspin>
inline decltype(auto) operator+(const std::array<T,nspin> &a, const std::array<T,nspin> &b){
  std::array<T,nspin> result;
#pragma unroll
  for(int i = 0; i < nspin; i++){
    result[i] = a[i].real() + b[i].real();
    result[i] = a[i].imag() + b[i].imag();
  }

  return result;
}

