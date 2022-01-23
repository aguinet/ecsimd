#ifndef ECSIMD_IFELSE_H
#define ECSIMD_IFELSE_H

#include <ecsimd/bignum.h>
#include <ecsimd/mgry.h>
#include <ecsimd/jacobian_curve_point.h>

#include <eve/function/if_else.hpp>

#include <cassert>

namespace ecsimd {

// For each bignum at idx i, return a[i] if mask[i] else b[i]
template <concepts::wide_bignum WBN>
auto if_else(cmp_res_t<WBN> mask, WBN& a, WBN& b)
{
  WBN ret;
  eve::detail::for_<0,1,bn_nlimbs<WBN>>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    eve::get<i>(ret) = eve::if_else(mask, eve::get<i>(a), eve::get<i>(b));
  });
  return ret;
}

template <concepts::wide_mgry_bignum WMBN>
auto if_else(cmp_res_t<typename WMBN::wide_bignum_type> mask, WMBN& a, WMBN& b)
{
  return WMBN{if_else(mask, a.wbn(), b.wbn())};
}

template <concepts::GFp GFP>
auto if_else(cmp_res_t<typename GFP::WBN> mask, GFP& a, GFP& b)
{
  return GFP{if_else(mask, a.wbn(), b.wbn())};
}

template <concepts::curve Curve>
auto if_else(
    cmp_res_t<curve_wide_bn_t<Curve>> mask,
    wide_jacobian_curve_point<Curve>& A,
    wide_jacobian_curve_point<Curve>& B)
{
  wide_jacobian_curve_point<Curve> Ret;
  Ret.x() = if_else(mask, A.x(), B.x());
  Ret.y() = if_else(mask, A.y(), B.y());
  Ret.z() = if_else(mask, A.z(), B.z());
  return Ret;
}

} // ecsimd

#endif
