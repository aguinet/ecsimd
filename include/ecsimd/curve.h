#ifndef ECSIMD_CURVE_H
#define ECSIMD_CURVE_H

#include <ecsimd/bignum.h>
#include <ecsimd/mgry.h>

namespace ecsimd {

namespace concepts {

// Weierstrass curves with a == -3
template <class T>
concept wst_curve_am3 = bignum<typename T::bn_type> &&
  bignum_cst<typename T::P>  && bignum_cst<typename T::B> &&
  bignum_cst<typename T::Gx> && bignum_cst<typename T::Gy>;

// Any possible curve
template <class T>
concept curve = wst_curve_am3<T>;

} // concepts

template <concepts::curve Curve>
using curve_bn_t = typename Curve::bn_type;

template <concepts::curve Curve>
using curve_wide_bn_t = wide_bignum<curve_bn_t<Curve>>;

template <concepts::curve Curve>
using curve_wide_mgry_bn_t = wide_mgry_bignum<curve_wide_bn_t<Curve>, typename Curve::P>;

} // ecsimd

#endif
