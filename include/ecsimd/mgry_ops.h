#ifndef ECSIMD_MGRY_OPS_H
#define ECSIMD_MGRY_OPS_H

#include <ecsimd/mgry.h>
#include <ecsimd/mgry_mul.h>
#include <ecsimd/modular.h>

namespace ecsimd {

template <concepts::wide_mgry_bignum WMBN>
WMBN mgry_add(WMBN const& a, WMBN const& b) {
  return WMBN{mod_add(a.wbn(), b.wbn(), WMBN::constants_type::wide_P)};
}

template <size_t Count, concepts::wide_mgry_bignum WMBN>
WMBN mgry_shift_left(WMBN const& a) {
  static_assert(Count > 0);
  WMBN ret(mod_shift_left_one(a.wbn(), WMBN::constants_type::wide_P));
  for (size_t i = 1; i < Count; ++i) {
    ret.wbn() = mod_shift_left_one(ret.wbn(), WMBN::constants_type::wide_P);
  }
  return ret;
}

template <concepts::wide_mgry_bignum WMBN>
WMBN mgry_sub(WMBN const& a, WMBN const& b) {
  return WMBN{mod_sub(a.wbn(), b.wbn(), WMBN::constants_type::wide_P)};
}

template <concepts::wide_mgry_bignum WMBN>
[[gnu::flatten]] WMBN mgry_mul(WMBN const& a, WMBN const& b) {
  auto m = mul(a.wbn(), b.wbn());
  return WMBN{details::mgry_reduce<typename WMBN::P_type>(m)};
  //return WMBN{details::mgry_mul<typename WMBN::P_type>(a.wbn(), b.wbn())};
}

template <concepts::wide_mgry_bignum WMBN>
[[gnu::flatten]] WMBN mgry_sqr(WMBN const& v) {
  auto s = square(v.wbn());
  return WMBN{details::mgry_reduce<typename WMBN::P_type>(s)};
  //return mgry_mul(v,v);
}

template <concepts::wide_bignum WBN, concepts::bignum_cst P>
wide_mgry_bignum<WBN, P> wide_mgry_bignum<WBN, P>::sqr() const {
  return mgry_sqr(*this);
}

// computes a**M*R [p]. Not safe from side channels leak of the exponent M.
template <concepts::wide_mgry_bignum WMBN>
WMBN mgry_pow(WMBN const& a, typename WMBN::bignum_type const& M) {
  using constants_t = typename WMBN::constants_type;
  using WBN = typename WMBN::wide_bignum_type;
  using limb_type = bn_limb_t<WBN>;
  constexpr size_t nbits = std::numeric_limits<limb_type>::digits;

  auto result = WMBN{WBN{constants_t::R_p}};
  const auto m = M.cbn();

  auto it_bzero = std::find_if(std::rbegin(m), std::rend(m), [](auto v) { return v != 0; }).base();
  if (it_bzero == std::begin(m)) {
    return result;
  }
  --it_bzero;

  auto base = a;
  for (auto itl = std::begin(m); itl != it_bzero; ++itl) {
    auto limb = *itl;
    for (size_t b = 0; b < nbits; ++b) {
      const auto lsb = limb & 1;
      limb >>= 1;
      if (lsb) {
        result = mgry_mul(result, base);
      }
      base = mgry_sqr(base);
    }
  }

  auto limb = *it_bzero;
  while (limb != 0) {
    const auto lsb = limb & 1;
    limb >>= 1;
    if (lsb) {
      result = mgry_mul(result, base);
    }
    if (limb == 0) break;
    base = mgry_sqr(base);
  }

  return result;
}

template <concepts::wide_mgry_bignum WMBN>
WMBN operator+(WMBN const& a, WMBN const& b) {
  return mgry_add(a,b);
}

template <concepts::wide_mgry_bignum WMBN>
WMBN operator-(WMBN const& a, WMBN const& b) {
  return mgry_sub(a,b);
}

template <concepts::wide_mgry_bignum WMBN>
WMBN operator*(WMBN const& a, WMBN const& b) {
  return mgry_mul(a,b);
}

} // ecsimd

#endif
