#ifndef ECSIMD_MGRY_CSTS_H
#define ECSIMD_MGRY_CSTS_H

#include <ecsimd/bignum.h>
#include <ctbignum/slicing.hpp>
#include <ctbignum/division.hpp>

namespace ecsimd {

template <concepts::wide_bignum WBN, concepts::bignum_cst P>
struct mgry_constants {
  using limb_type = bn_limb_t<WBN>;
  using BN = typename WBN::value_type;
  static constexpr size_t nlimbs = bn_nlimbs<WBN>;
  static constexpr auto R   = cbn::detail::unary_encoding<nlimbs, nlimbs + 1, limb_type>();
  static constexpr auto Rsq = cbn::detail::unary_encoding<2 * nlimbs, 2 * nlimbs + 1, limb_type>();

  static constexpr auto p = P::value.cbn();

  static constexpr auto R_p   = BN::from(cbn::div(R, p).remainder);
  static constexpr auto Rsq_p = BN::from(cbn::div(Rsq, p).remainder);

  using cbn_type = typename BN::cbn_type;
  static constexpr auto Pm1_by_R_p = BN::from(cbn::div((p-cbn_type{1})*R, p).remainder);

  static const WBN wide_P;
  static const wbn_zext_t<WBN> wide_P_zext;
};


template <concepts::wide_bignum WBN, concepts::bignum_cst P>
const WBN mgry_constants<WBN, P>::wide_P = WBN{P::value};

template <concepts::wide_bignum WBN, concepts::bignum_cst P>
const wbn_zext_t<WBN> mgry_constants<WBN, P>::wide_P_zext = zext_u32x64(WBN{P::value});

} // ecsimd

#endif
