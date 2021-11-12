#ifndef ECSIMD_GFP_H
#define ECSIMD_GFP_H

#include <ecsimd/bignum.h>
#include <ecsimd/mgry.h>
#include <ecsimd/mgry_ops.h>

#include <eve/function/any.hpp>

#include <ctbignum/addition.hpp>
#include <ctbignum/division.hpp>

#include <optional>

namespace ecsimd {

template <concepts::bignum_cst P>
struct GFp
{
  using BN = bn_t<P>;
  using WBN = wide_bignum<BN>;
  using WMBN = wide_mgry_bignum<WBN, P>;

  static auto mgry_inverse(WMBN const& v) {
    return mgry_pow(v, P_pow_m2);
  }

  static std::optional<WMBN> mgry_sqrt(WMBN const& v) {
    const auto ret = mgry_pow(v, P_pow_sqrt);
    // TODO: verify it exists using the Jacobi symbol
    // (https://en.wikipedia.org/wiki/Jacobi_symbol#Calculating_the_Jacobi_symbol)
    if (eve::any(mgry_sqr(ret).wbn() != v.wbn())) {
      return {};
    }
    return {ret};
  }

  static WMBN mgry_opposite(WMBN const& v) {
    using constants_type = typename WMBN::constants_type;
    WMBN ret = mgry_sub(v, WMBN::R());
    return mgry_sub(WMBN{WBN{constants_type::Pm1_by_R_p}}, ret);
  }

private:
  using cbn_type = typename BN::cbn_type;

  static constexpr auto P_cbn = P::value.cbn();
  static constexpr auto P_pow_m2 = BN::from(cbn::subtract_ignore_carry(
      P_cbn,cbn_type{2}));
  // See https://www.rieselprime.de/ziki/Modular_square_root
  // P % 4 == 3
  static_assert((P_cbn[0] & 3) == 3);
  static constexpr auto P_pow_sqrt = BN::from(
    cbn::detail::first<bn_nlimbs<BN>>(
      cbn::shift_right(P_cbn+cbn_type{1},2)));
};

template <concepts::wide_bignum WBN, concepts::bignum_cst P>
wide_mgry_bignum<WBN, P> wide_mgry_bignum<WBN, P>::inverse() const {
  return GFp<P>::mgry_inverse(*this);
}

template <concepts::wide_bignum WBN, concepts::bignum_cst P>
std::optional<wide_mgry_bignum<WBN, P>> wide_mgry_bignum<WBN, P>::sqrt() const {
  return GFp<P>::mgry_sqrt(*this);
}

} // ecsimd

#endif
