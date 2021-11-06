#ifndef FPSIMD_BIGINT_SUB_H
#define FPSIMD_BIGINT_SUB_H

#include <eve/wide.hpp>
#include <ecsimd/bignum.h>

#include <tuple>

namespace ecsimd {

template <concepts::wide_bignum WBN>
auto sub(WBN const& a, WBN const& b)
{
  using limb_type = bn_limb_t<WBN>;
  using C = typename WBN::cardinal_type;
  constexpr auto nlimbs = bn_nlimbs<WBN>;

  cmp_res_t<WBN> carry_mask;
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
auto sub_no_carry(WBN const& a, WBN const& b) {
  return std::get<0>(sub(a,b));
}

// return a-p if a > p, else a
template <size_t nlimbs_out, concepts::wide_bignum WBN, concepts::cmp_res<WBN>... Carry>
auto sub_if_above(WBN const& a, WBN const& p, Carry... additional_carries)
{
  using limb_type = bn_limb_t<WBN>;
  using cardinal = eve::cardinal_t<WBN>;
  constexpr auto nlimbs = bn_nlimbs<WBN>;

  WBN asub;
  cmp_res_t<WBN> carry_mask;
  std::tie(asub, carry_mask) = sub(a, p);
  if constexpr (sizeof...(additional_carries) > 0) {
    carry_mask = (carry_mask || (additional_carries || ...));
  }

  eve::wide<bignum<limb_type, nlimbs_out>, cardinal> ret;
  eve::detail::for_<0,1,nlimbs_out>([&](auto i_) {
    constexpr auto i = decltype(i_)::value;
    const auto va = eve::get<i>(a);
    const auto vasub = eve::get<i>(asub);
    const auto mask = eve::broadcast(carry_mask, eve::index<i>);

    eve::get<i>(ret) = eve::if_else(mask, va, vasub);
  });
  return ret;
}

template <concepts::wide_bignum WBN, concepts::cmp_res<WBN>... Carry>
auto sub_if_above(WBN const& a, WBN const& p, Carry... additional_carries) {
  return sub_if_above<bn_nlimbs<WBN>, WBN>(a,p,
      std::forward<Carry>(additional_carries)...);
}

} // ecsimd

#endif
