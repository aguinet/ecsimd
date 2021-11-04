#ifndef ECSIMD_CMP_H
#define ECSIMD_CMP_H

#include <eve/wide.hpp>
#include <ecsimd/bignum.h>
#include <ecsimd/sub.h>

namespace ecsimd {

template <concepts::wide_bignum WBN>
using cmp_ret_t = eve::wide<bn_limb_t<WBN>, eve::cardinal_t<WBN>>;

template <concepts::wide_bignum WBN>
static cmp_ret_t<WBN> cmp_lt(WBN const& a, WBN const& b) {
  return std::get<1>(sub(a,b));
}

template <concepts::wide_bignum WBN>
static cmp_ret_t<WBN> cmp_gt(WBN const& a, WBN const& b) {
  return cmp_lt(b,a);
}

template <concepts::wide_bignum WBN>
static cmp_ret_t<WBN> cmp_lte(WBN const& a, WBN const& b) {
  return ~cmp_gt(a,b);
}

template <concepts::wide_bignum WBN>
static cmp_ret_t<WBN> cmp_gte(WBN const& a, WBN const& b) {
  return ~cmp_lt(a,b);
}

} // ecsimd

template <ecsimd::concepts::wide_bignum WBN>
static ecsimd::cmp_ret_t<WBN> operator<(WBN const& a, WBN const& b) {
  return ecsimd::cmp_lt(a,b);
}

template <ecsimd::concepts::wide_bignum WBN>
static ecsimd::cmp_ret_t<WBN> operator>(WBN const& a, WBN const& b) {
  return ecsimd::cmp_gt(a,b);
}

template <ecsimd::concepts::wide_bignum WBN>
static ecsimd::cmp_ret_t<WBN> operator<=(WBN const& a, WBN const& b) {
  return ecsimd::cmp_lte(a,b);
}

template <ecsimd::concepts::wide_bignum WBN>
static ecsimd::cmp_ret_t<WBN> operator>=(WBN const& a, WBN const& b) {
  return ecsimd::cmp_gte(a,b);
}

#endif
