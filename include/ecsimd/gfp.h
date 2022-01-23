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

template <concepts::wide_bignum WBN_, concepts::bignum_cst P>
struct GFp
{
  using P_type = P;
  using WBN = WBN_;
  using BN = typename WBN::value_type;
  using WMBN = wide_mgry_bignum<WBN, P>;

  GFp() = default;
  GFp(WMBN const& n):
    n_(n)
  { }

  static GFp one() {
    return GFp{WMBN::R()};
  }

  static GFp from_classical(WBN const& n) {
    return {WMBN::from_classical(n)};
  }

  WBN to_classical() const {
    return n_.to_classical();
  }

  GFp inverse() const {
    return {mgry_pow(n_, P_pow_m2)};
  }

  std::optional<GFp> sqrt() const {
    const auto ret = GFp{mgry_pow(n_, P_pow_sqrt)};
    // TODO: verify it exists using the Jacobi symbol
    // (https://en.wikipedia.org/wiki/Jacobi_symbol#Calculating_the_Jacobi_symbol)
    if (eve::any(ret.sqr().wbn() != wbn())) {
      return {};
    }
    return {ret};
  }

  GFp sqr() const {
    return {mgry_sqr(n_)};
  }

  GFp opposite() const {
    using constants_type = typename WMBN::constants_type;
    WMBN ret = mgry_sub(n_, WMBN::R());
    return GFp{mgry_sub(WMBN{WBN{constants_type::Pm1_by_R_p}}, ret)};
  }

  auto operator<=>(GFp const& o) const {
    return n_ <=> o.n_;
  }

  auto const& wbn() const { return n_.wbn(); }
  auto& wbn() { return n_.wbn(); }

  auto const& wmbn() const { return n_; }
  auto& wmbn() { return n_; }

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

  WMBN n_;
};

namespace concepts {
template <class T>
concept GFp = std::same_as<T, GFp<typename T::WBN, typename T::P_type>>;
} // concepts

template <concepts::GFp GFP>
GFP operator+(GFP const& a, GFP const& b) {
  return GFP{mgry_add(a.wmbn(),b.wmbn())};
}

template <concepts::GFp GFP>
GFP operator-(GFP const& a, GFP const& b) {
  return GFP{mgry_sub(a.wmbn(),b.wmbn())};
}

template <concepts::GFp GFP>
GFP operator*(GFP const& a, GFP const& b) {
  return GFP{mgry_mul(a.wmbn(),b.wmbn())};
}

template <size_t Count, concepts::GFp GFP>
GFP gfp_shift_left(GFP const& a) {
  return GFP{mgry_shift_left<Count>(a.wmbn())};
}

} // ecsimd

#endif
