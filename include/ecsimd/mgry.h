#ifndef ECSIMD_MGRY_H
#define ECSIMD_MGRY_H

#include <ecsimd/bignum.h>
#include <ecsimd/mgry_csts.h>
#include <ecsimd/shift.h>
#include <eve/product_type.hpp>

namespace ecsimd {

namespace details {
template <concepts::bignum_cst P, concepts::wide_bignum WBN>
WBN mgry_mul(WBN const& a, WBN const& b);
template <concepts::bignum_cst P_type, concepts::wide_bignum WBN>
auto mgry_reduce(WBN const& a);
} // details

template <concepts::bignum_cst P>
constexpr auto to_mgry(bn_t<P> const& v) {
  using BN = bn_t<P>;
  using WBN = wide_bignum<BN>;

  const auto v_cbn = v.cbn();
  const auto p = P::value.cbn();
  return BN::from(cbn::div(v_cbn * mgry_constants<WBN, P>::R, p).remainder);
}

template <concepts::wide_bignum WBN, concepts::bignum_cst P>
struct wide_mgry_bignum
{
  using wide_bignum_type = WBN;
  using bignum_type = typename wide_bignum_type::value_type;
  using P_type = P;

  using constants_type = mgry_constants<wide_bignum_type, P_type>;

  wide_mgry_bignum() = default;

  wide_mgry_bignum(WBN const& n):
    n_(n)
  { }

  static wide_mgry_bignum R() {
    return wide_mgry_bignum{WBN{constants_type::R_p}};
  }

  static wide_mgry_bignum from_classical(WBN const& n) {
    const auto m = mul(n, WBN{constants_type::Rsq_p});
    return wide_mgry_bignum{details::mgry_reduce<P_type>(m)};
  }

  wide_mgry_bignum sqr() const;
  wide_mgry_bignum inverse() const;
  std::optional<wide_mgry_bignum> sqrt() const;

  WBN to_classical() const {
    auto zext = pad<bn_nlimbs<bignum_type>>(n_);
    return details::mgry_reduce<P_type>(zext);
  }

  WBN const& wbn() const { return n_; }
  WBN& wbn() { return n_; }

  auto operator<=>(wide_mgry_bignum const& o) const {
    return wbn() <=> o.wbn();
  }

private:
  WBN n_;
};

namespace concepts {
template <class T>
concept wide_mgry_bignum = std::same_as<T, wide_mgry_bignum<typename T::wide_bignum_type, typename T::P_type>>;
} // concepts

} // ecsimd

#endif
