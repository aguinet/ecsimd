#ifndef FPSIMD_BIGINT_ADD_H
#define FPSIMD_BIGINT_ADD_H

#include <eve/wide.hpp>
#include <ecsimd/bignum.h>
#include <ecsimd/compat.h>

#include <tuple>

namespace ecsimd {

template <concepts::wide_bignum WBN>
static auto add(WBN const& a, WBN const& b)
{
  constexpr auto nlimbs = bn_nlimbs<WBN>;

  cmp_res_t<WBN> carry_mask;
  WBN ret;
  eve::detail::for_<0, 1, nlimbs>([&](auto i_) ECSIMD_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    const auto aa = eve::get<i>(a);
    const auto bb = eve::get<i>(b);
    const auto sum = aa + bb;
    if constexpr (i == 0) {
      carry_mask = sum < aa;
      eve::get<i>(ret) = sum;
    }
    else {
      const auto res = sum - carry_mask.mask();
      carry_mask = (sum < aa) || (res < sum);
      eve::get<i>(ret) = res;
    }
  });

  return std::make_tuple(ret, carry_mask);
}

template <concepts::wide_bignum WBN>
static auto add_no_carry(WBN const& a, WBN const& b)
{
  return std::get<0>(add(a,b));
}

} // ecsimd

#endif
