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
  eve::detail::for_<0,1,bn_nlimbs<WBN>>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
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

template <size_t N, concepts::wide_bignum WBN>
const auto pad(WBN const& v) {
  using limb_type = bn_limb_t<WBN>;
  constexpr auto v_nlimbs = bn_nlimbs<WBN>;
  constexpr auto ret_nlimbs = v_nlimbs+N;
  using ret_ty = wide_bignum<bignum<limb_type, ret_nlimbs>>;

  ret_ty ret;
  eve::detail::for_<0,1,v_nlimbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    eve::get<i>(ret) = eve::get<i>(v);
  });
  eve::detail::for_<v_nlimbs,1,ret_nlimbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    eve::get<i>(ret) = eve::zero(eve::as(eve::get<i>(ret)));
  });
  return ret;
}

template <size_t RetLimbs, size_t ShiftBy, concepts::wide_bignum WBN>
auto limb_shift_left(WBN const& v)
{
  using limb_type = bn_limb_t<WBN>;
  constexpr auto v_nlimbs = bn_nlimbs<WBN>;
  using ret_ty = wide_bignum<bignum<limb_type, RetLimbs>>;

  if constexpr (ShiftBy >= RetLimbs || ShiftBy >= v_nlimbs) {
    return eve::zero(eve::as<ret_ty>());
  }

  ret_ty ret;
  eve::detail::for_<0,1,ShiftBy>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    eve::get<i>(ret) = eve::zero(eve::as(eve::get<i>(ret)));
  });
  eve::detail::for_<ShiftBy,1,ShiftBy+v_nlimbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    eve::get<i>(ret) = eve::get<i-ShiftBy>(v);
  });
  eve::detail::for_<ShiftBy+v_nlimbs, 1, RetLimbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    eve::get<i>(ret) = eve::zero(eve::as(eve::get<i>(ret)));
  });
  return ret;
}

template <size_t ShiftBy, concepts::wide_bignum WBN>
auto limb_shift_right(WBN const& v)
{
  using limb_type = bn_limb_t<WBN>;
  constexpr auto v_nlimbs = bn_nlimbs<WBN>;
  static_assert(ShiftBy < v_nlimbs);
  constexpr auto ret_limbs = v_nlimbs-ShiftBy;
  using ret_ty = wide_bignum<bignum<limb_type, ret_limbs>>;


  ret_ty ret;
  eve::detail::for_<0,1,ret_limbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    eve::get<i>(ret) = eve::get<i+ShiftBy>(v);
  });
  return ret;
}

} // ecsimd

#endif
