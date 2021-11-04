#ifndef FPSIMD_BIGINT_SUB_H
#define FPSIMD_BIGINT_SUB_H

#include <eve/wide.hpp>
#include <ecsimd/bignum.h>

#include <tuple>

namespace ecsimd {

template <concepts::wide_bignum WBN>
static auto sub(WBN const& a, WBN const& b)
{
  using limb_type = bn_limb_t<WBN>;
  using C = typename WBN::cardinal_type;
  constexpr auto nlimbs = bn_nlimbs<WBN>;

  eve::logical<eve::wide<limb_type, C>> carry_mask;
  WBN ret;
  eve::detail::for_<0, 1, nlimbs>([&](auto i_) {
    constexpr auto i = decltype(i_)::value;
    const auto aa = eve::get<i>(a);
    const auto bb = eve::get<i>(b);
    const auto diff = aa - bb;
    if constexpr (i == 0) {
      carry_mask = diff > aa;
      eve::get<i>(ret) = diff;
    }
    else {
      const auto res = diff + carry_mask.mask();
      carry_mask = (diff > aa) || (res > diff);
      eve::get<i>(ret) = res;
    }
  });

  return std::make_tuple(ret, carry_mask);
}

template <concepts::wide_bignum WBN>
static auto sub_no_carry(WBN const& a, WBN const& b) {
  return std::get<0>(sub(a,b));
}

} // ecsimd

#endif
