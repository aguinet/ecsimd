#ifndef ECSIMD_SELECT_H
#define ECSIMD_SELECT_H

#include <ecsimd/bignum.h>
#include <ecsimd/mgry.h>
#include <ecsimd/jacobian_curve_point.h>

#include <eve/function/swap_if.hpp>

#include <cassert>

namespace ecsimd {

// For each bignum at idx i, swap (a[i],b[i]) if mask[i] is true
template <concepts::wide_bignum WBN>
void swap_if(cmp_res_t<WBN> mask, WBN& a, WBN& b)
{
  eve::detail::for_<0,1,bn_nlimbs<WBN>>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    eve::swap_if(mask, eve::get<i>(a), eve::get<i>(b));
  });
}

template <concepts::wide_mgry_bignum WMBN>
void swap_if(cmp_res_t<typename WMBN::wide_bignum_type> mask, WMBN& a, WMBN& b)
{
  swap_if(mask, a.wbn(), b.wbn());
}

template <concepts::curve Curve>
void swap_if(
    cmp_res_t<curve_wide_bn_t<Curve>> mask,
    wide_jacobian_curve_point<Curve>& A,
    wide_jacobian_curve_point<Curve>& B)
{
  swap_if(mask, A.x(), B.x());
  swap_if(mask, A.y(), B.y());
  swap_if(mask, A.z(), B.z());
}

template <concepts::curve Curve>
void swap_if_same_z(
    cmp_res_t<curve_wide_bn_t<Curve>> mask,
    wide_jacobian_curve_point<Curve>& A,
    wide_jacobian_curve_point<Curve>& B)
{
  assert(eve::all(A.z().wbn() == B.z().wbn()));
  swap_if(mask, A.x(), B.x());
  swap_if(mask, A.y(), B.y());
}

} // ecsimd

#endif
