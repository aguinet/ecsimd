#ifndef ECSIMD_MODULAR_H
#define ECSIMD_MODULAR_H

#include <ecsimd/add.h>
#include <ecsimd/sub.h>
#include <ecsimd/shift.h>

namespace ecsimd {

template <concepts::wide_bignum WBN>
auto mod_add(WBN const& a, WBN const& b, WBN const& p)
{
  const auto [sum, carry_add] = add(a,b);
  return sub_if_above(sum, p, !carry_add);
}

template <concepts::wide_bignum WBN>
auto mod_shift_left_one(WBN const& a, WBN const& p)
{
  const auto [shifted, carry] = shift_left_one(a);
  return sub_if_above(shifted, p, !carry);
}

template <concepts::wide_bignum WBN>
auto mod_sub(WBN const& a, WBN const& b, WBN const& p)
{
  WBN diff;
  cmp_res_t<WBN> carry_sub;
  std::tie(diff, carry_sub) = sub(a,b);
  const WBN diff_add = add_no_carry(diff, p);

  WBN ret;
  eve::detail::for_<0,1,bn_nlimbs<WBN>>([&](auto i_) ECSIMD_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    const auto vd  = eve::get<i>(diff);
    const auto vda = eve::get<i>(diff_add);

    eve::get<i>(ret) = eve::if_else(carry_sub, vda, vd);
  });
  return ret;
}

} // ecsimd

#endif
