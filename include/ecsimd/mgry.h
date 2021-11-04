#ifndef ECSIMD_MGRY_H
#define ECSIMD_MGRY_H

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
};

namespace details {
/*
template <concepts::bignum_cst P, concepts::wide_bignum WBN>
WBN mgry_mul(WBN const& a, WBN const& b);*/

template <concepts::bignum_cst P, concepts::wide_bignum WBN>
WBN mgry_reduce(WBN const& a);
}

template <concepts::wide_bignum WBN, concepts::bignum_cst P>
struct wide_mgry_bignum
{
  using wide_bignum_type = WBN;
  using P_type = P;

  using constants_type = mgry_constants<wide_bignum_type, P_type>;

  wide_mgry_bignum(WBN const& n):
    n_(n)
  { }

  /*
  static wide_mgry_bignum from_classical(WBN const& n) {
    return wide_mgry_bignum{details::mgry_mul<P_type>(n, WBN{constants_type::Rsq_p})};
  }*/

  WBN to_classical() const {
    return details::mgry_reduce<P_type>(n_);
  }

  WBN const& wbn() const { return n_; }

private:
  WBN n_;
};

namespace concepts {
template <class T>
concept wide_mgry_bignum = std::same_as<T, wide_mgry_bignum<typename T::wide_bignum_type, typename T::P_type>>;
} // concepts

// Operations

/*
template <concepts::wide_mgry_bignum WMBN>
auto mul(WMBN const& a, WMBN const& b) {
  return WMBN{details::mgry_mul<typename WMBN::P_type>(a.bn(), b.bn())};
}*/

} // ecsimd::mgry

#endif