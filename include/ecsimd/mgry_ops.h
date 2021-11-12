#ifndef ECSIMD_MGRY_OPS_H
#define ECSIMD_MGRY_OPS_H

#include <ecsimd/mgry.h>
#include <ecsimd/mgry_mul.h>
#include <ecsimd/modular.h>

namespace ecsimd {

template <concepts::wide_mgry_bignum WMBN>
WMBN mgry_add(WMBN const& a, WMBN const& b) {
  return WMBN{mod_add(a.wbn(), b.wbn(), WMBN::wide_P)};
}

template <concepts::wide_mgry_bignum WMBN>
WMBN mgry_sub(WMBN const& a, WMBN const& b) {
  return WMBN{mod_sub(a.wbn(), b.wbn(), WMBN::wide_P)};
}

template <concepts::wide_mgry_bignum WMBN>
WMBN mgry_mul(WMBN const& a, WMBN const& b) {
  return WMBN{details::mgry_mul<WMBN::P_type>(a.wbn(), b.wbn())};
}

template <concepts::wide_mgry_bignum WMBN>
WMBN mgry_sqr(WMBN const& v) {
  // TODO: optimize this
  return mgry_mul(v,v);
}

// computes a**m*R [p]
template <concepts::wide_mgry_bignum WMBN>
WMBN mgry_pow(WMBN const& a, typename WMBN::bignum_type const& M) {
  using constants_t = typename WMBN::constants_type;
  using WBN = typename WMBN::wide_bignum_type;

  auto result = WMBN{WBN{constants_t::R_p}};
  const auto m = M.cbn();
  auto base = a;
  for (auto const limb: m) {
    while (limb != 0) {
      const auto lsb = limb & 1;
      limb >>= 1;
      if (lsb) {
        result = mgry_mul(result, base);
        if (limb == 0) break;
      }
      base = mgry_sqr(base);
    }
  }
  return result;
}

} // ecsimd

#endif
