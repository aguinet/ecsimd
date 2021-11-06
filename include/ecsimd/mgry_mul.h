#ifndef ECSIMD_MGRY_MUL_H
#define ECSIMD_MGRY_MUL_H

#include <ecsimd/mul.h>
#include <ecsimd/sub.h>
#include <ecsimd/utility.h>
#include <ecsimd/cmp.h>

#include <ctbignum/gcd.hpp>
#include <ctbignum/io.hpp>

#include <eve/wide.hpp>
#include <eve/detail/meta.hpp>
#include <eve/traits/cardinal.hpp>

#include <utility>
#include <limits>

namespace ecsimd::details {

template <concepts::bignum_cst P_type>
class mgry_constants {
  static constexpr auto P = P_type::value;
  static constexpr auto P_cbn = P.cbn();
  using BN = std::decay_t<decltype(P)>;
  using limb_type = bn_limb_t<BN>;

  static constexpr auto inv = cbn::mod_inv(
    array_to_integer_sequence_t<P_cbn>{},
    std::integer_sequence<limb_type, 0, 1>{}); // 2**32

public:
  static constexpr limb_type mprime = -inv[0];
};

template <concepts::bignum_cst P_type, concepts::wide_bignum WBN>
__attribute__((noinline)) WBN mgry_mul(WBN const& x, WBN const& y) {
  using limb_type = bn_limb_t<WBN>;
  using dbl_limb_type = eve::detail::upgrade_t<limb_type>;
  constexpr auto limb_bits = std::numeric_limits<limb_type>::digits;
  constexpr auto nlimbs = bn_nlimbs<WBN>;
  using cardinal = eve::cardinal_t<WBN>;
  using bignum_p1 = bignum<limb_type, nlimbs+1>;
  using wide_bignum_p1 = wide_bignum<bignum_p1>;

  constexpr auto P = P_type::value;
  constexpr auto PNp1 = bignum_p1::from(cbn::detail::pad<1>(P.cbn()));
  constexpr auto mprime = mgry_constants<P_type>::mprime;

  using wide_limb_type = eve::wide<limb_type, cardinal>;
  using wide_dbl_limb_type = eve::wide<dbl_limb_type, cardinal>;

  const wide_limb_type wide_mprime(mprime);

  const auto wide_dbl_low_mask = wide_dbl_limb_type{static_cast<dbl_limb_type>(
          std::numeric_limits<limb_type>::max())};

  auto A = eve::zero(eve::as<wide_bignum_p1>());
  eve::detail::for_<0, 1, nlimbs>([&](auto i_) {
      constexpr auto i = decltype(i_)::value;

      const auto u_i = (eve::get<0>(A) + eve::get<i>(x) * eve::get<0>(y)) * wide_mprime;

      auto k = eve::zero(eve::as<wide_dbl_limb_type>());
      auto k2 = eve::zero(eve::as<wide_dbl_limb_type>());

      auto z = mul_wide(eve::get<0>(y), eve::get<i>(x));
      z += eve::convert(eve::get<0>(A), eve::as<dbl_limb_type>());
      z += k;

      auto z2 = mul_wide(wide_limb_type{kumi::get<0>(P)}, u_i);
      z2 += z & wide_dbl_low_mask;
      z2 += k2;

      k = z >> limb_bits;
      k2 = z2 >> limb_bits;

      eve::detail::for_<1,1,nlimbs>([&](auto j_) {
        constexpr auto j = decltype(j_)::value;

        auto t = mul_wide(eve::get<j>(y), eve::get<i>(x));
        t += eve::convert(eve::get<j>(A), eve::as<dbl_limb_type>());
        t += k;

        auto t2 = mul_wide(wide_limb_type{kumi::get<j>(P)}, u_i);
        t2 += t & wide_dbl_low_mask;
        t2 += k2;
        eve::get<j-1>(A) = eve::convert(t2, eve::as<limb_type>());

        k = t >> limb_bits;
        k2 = t2 >> limb_bits;
      });

      const auto tmp = eve::convert(eve::get<nlimbs>(A), eve::as<dbl_limb_type>()) + k + k2;
      // TODO: these two lines could be optimized with shuffles + stores
      eve::get<nlimbs-1>(A) = eve::convert(tmp, eve::as<limb_type>());
      eve::get<nlimbs>(A) = eve::convert(tmp >> limb_bits, eve::as<limb_type>());
  });

  return sub_if_above<nlimbs>(A, wide_bignum_p1{PNp1});
}

} // ecsimd::details

#endif
