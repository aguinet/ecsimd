#ifndef ECSIMD_MGRY_MUL_H
#define ECSIMD_MGRY_MUL_H

#include <ecsimd/mgry_csts.h>
#include <ecsimd/add.h>
#include <ecsimd/mul.h>
#include <ecsimd/sub.h>
#include <ecsimd/shift.h>
#include <ecsimd/utility.h>
#include <ecsimd/cmp.h>

#include <ctbignum/gcd.hpp>
#include <ctbignum/io.hpp>

#include <eve/wide.hpp>
#include <eve/detail/meta.hpp>
#include <eve/traits/cardinal.hpp>

#include <utility>
#include <limits>
#include <iostream>

namespace ecsimd::details {

template <concepts::bignum_cst P_type, typename Cardinal>
class mgry_mul_constants {
  static constexpr auto P = P_type::value;
  static constexpr auto P_cbn = P.cbn();
  using BN = std::decay_t<decltype(P)>;
  using limb_type = bn_limb_t<BN>;
  static constexpr auto nlimbs = bn_nlimbs<BN>;

  static constexpr auto inv = cbn::mod_inv(
    array_to_integer_sequence_t<P_cbn>{},
    std::integer_sequence<limb_type, 0, 1>{}); // 2**32

  using bignum_p1 = bignum<limb_type, nlimbs+1>;

public:
  static constexpr limb_type mprime = -inv[0];
  static constexpr auto PNp1 = bignum_p1::from(cbn::detail::pad<1>(P.cbn()));
  static const eve::wide<bn_limb_t<P_type>, Cardinal> wide_mprime;
  static const eve::wide<eve::detail::upgrade_t<bn_limb_t<P_type>>, Cardinal> wide_mprime_zext;
};

template <concepts::bignum_cst P_type, typename Cardinal>
const eve::wide<bn_limb_t<P_type>, Cardinal> mgry_mul_constants<P_type, Cardinal>::wide_mprime = eve::wide<bn_limb_t<P_type>, Cardinal>{mgry_mul_constants<P_type, Cardinal>::mprime};

template <concepts::bignum_cst P_type, typename Cardinal>
const eve::wide<eve::detail::upgrade_t<bn_limb_t<P_type>>, Cardinal> mgry_mul_constants<P_type, Cardinal>::wide_mprime_zext = eve::convert(eve::wide<bn_limb_t<P_type>, Cardinal>{mgry_mul_constants<P_type, Cardinal>::mprime}, eve::as<eve::detail::upgrade_t<bn_limb_t<P_type>>>());

template <concepts::wide_bignum WBN>
static auto add_no_carry_u32_zext(WBN const& a, WBN const& b)
{
  constexpr auto nlimbs = bn_nlimbs<WBN>;
  using limb_type = bn_limb_t<WBN>;
  using half_limb_type = eve::detail::downgrade_t<limb_type>;
  using cardinal = eve::cardinal_t<WBN>;
  constexpr auto half_nbits = std::numeric_limits<half_limb_type>::digits;

  WBN ret;
  eve::detail::for_<0, 1, nlimbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    const auto aa = eve::get<i>(a);
    const auto bb = eve::get<i>(b);
    const auto sum = aa + bb;
    eve::get<i>(ret) = sum;
  });

  // Normalize to 32 bits limbs zexted to 64 bits
  const eve::wide<limb_type, cardinal> low_mask(std::numeric_limits<half_limb_type>::max());
  eve::detail::for_<1, 1, nlimbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    const auto prev = eve::get<i-1>(ret);
    eve::get<i>(ret) += prev >> half_nbits;
    eve::get<i-1>(ret) &= low_mask;
  });
  // TODO: is this necessary?
  eve::get<nlimbs-1>(ret) &= low_mask;

  return ret;
}

template <concepts::bignum_cst P_type, concepts::wide_bignum WBN>
__attribute__((flatten)) auto mgry_reduce(WBN const& a)
{
  static_assert(bn_nlimbs<WBN> == 2*bn_nlimbs<P_type>);
  static_assert(std::is_same_v<bn_limb_t<WBN>, bn_limb_t<P_type>>);

  using limb_type = bn_limb_t<WBN>;
  using half_limb_type = eve::detail::downgrade_t<limb_type>;
  using cardinal = eve::cardinal_t<WBN>;
  using half_P_type = remap_limb_t<P_type, half_limb_type>;
  using mul_csts = mgry_mul_constants<half_P_type, cardinal>;
  constexpr auto nlimbs = bn_nlimbs<P_type>;
  constexpr auto half_nlimbs = bn_nlimbs<half_P_type>;

  using wide_half_type = wide_bignum<bignum<half_limb_type, half_nlimbs>>;
  using wide_type = wide_bignum<bignum<limb_type, nlimbs>>;

  const auto a_zext = zext_u32x64(a);
  auto accum = pad<1>(a_zext);
  constexpr auto nlimbs_accum = bn_nlimbs<decltype(accum)>;

  // TODO: make sure these converts happens at compile time
  const auto& wide_P_zext = mgry_constants<wide_type, P_type>::wide_P_zext;
  const auto& wide_mprime = mul_csts::wide_mprime_zext;
  const eve::wide<limb_type, cardinal> low_mask(std::numeric_limits<half_limb_type>::max());

  eve::detail::for_<0,1,half_nlimbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    // '& low_mask' could be set to mullow, but limb_mul_zext already only considers the lowest 32bits
    auto prod = limb_mul_zext(wide_P_zext, mullow(eve::get<i>(accum), wide_mprime));
    auto prod2 = limb_shift_left<nlimbs_accum, i>(prod);
    accum = add_no_carry_u32_zext(accum, prod2);
  });

  auto result = trunc_u64x32(pad<1>(limb_shift_right<half_nlimbs>(accum)));
  const auto wide_P_pad1 = pad<1>(mgry_constants<wide_type, P_type>::wide_P);
  return sub_if_above<bn_nlimbs<P_type>>(result, wide_P_pad1);
}

} // ecsimd::details

#endif
