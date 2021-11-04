#ifndef ECSIMD_TESTS_H
#define ECSIMD_TESTS_H

#include <ecsimd/bignum.h>
#include <eve/wide.hpp>
#include <array>

namespace ecsimd {

template <concepts::wide_bignum WBN>
static auto wide_bignum_set1(std::array<uint8_t, sizeof(typename WBN::value_type)> const& bytes) {
  const auto BN = bn_from_bytes_BE<typename WBN::value_type>(bytes);
  return WBN{BN};
}

} // ecsimd

#endif
