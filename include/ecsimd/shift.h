#ifndef ECSIMD_SHIFT_H
#define ECSIMD_SHIFT_H

#include <ecsimd/bignum.h>
#include <ecsimd/utility.h>

#include <concepts>
#include <tuple>
#include <bit>

namespace ecsimd {

template <concepts::wide_bignum WBN>
auto shift_left_one(WBN const& a) {
  cmp_res_t<WBN> carry;
  WBN ret;
  eve::detail::for_<0,1,bn_nlimbs<WBN>>([&](auto i_) {
    constexpr auto i = decltype(i_)::value;

    const auto l = eve::get<i>(a);
    auto shifted = l << 1;
    if constexpr (i > 0) {
      shifted |= carry.mask() & 1U;
    }
    carry = wide_uasr(l, 31);
    eve::get<i>(ret) = shifted;
  });
  return std::make_tuple(ret, carry);
}

} // ecsimd

#endif
