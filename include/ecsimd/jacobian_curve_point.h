#ifndef ECSIMD_JACOBIAN_CURVE_POINT_H
#define ECSIMD_JACOBIAN_CURVE_POINT_H

#include <ecsimd/curve.h>
#include <ecsimd/curve_point.h>
#include <ecsimd/gfp.h>

namespace ecsimd {

// TODO: type erasure, only keep the bignum type as a template argument
template <concepts::curve Curve>
struct wide_jacobian_curve_point {
  using curve_type = Curve;
  using bignum_type = curve_bn_t<Curve>;
  using WBN  = curve_wide_bn_t<Curve>;
  using WMBN = curve_wide_mgry_bn_t<Curve>;
  using wide_curve_point_t = wide_curve_point<Curve>;
  using gfp = GFp<typename Curve::P>;

  wide_jacobian_curve_point() = default;
  wide_jacobian_curve_point(wide_jacobian_curve_point const&) = default;

  // affine_x = x / z**2
  // affine_y = y / z**3
  // We also use the montgomery representation, so z == R
  static wide_jacobian_curve_point from_affine(wide_curve_point_t const& pt) {
    wide_jacobian_curve_point ret;
    ret.x_ = WMBN::from_classical(pt.x());
    ret.y_ = WMBN::from_classical(pt.y());
    ret.z_ = WMBN::R();
    return ret;
  }

  wide_curve_point_t to_affine() const {
    const auto invZ = z_.inverse();
    const auto invZ2 = invZ.sqr();
    const auto invZ3 = invZ2 * invZ;

    wide_curve_point_t ret;
    ret.x() = (x_ * invZ2).to_classical();
    ret.y() = (y_ * invZ3).to_classical();
    return ret;
  }

  auto operator==(wide_jacobian_curve_point const& o) const {
    return x().wbn() == o.x().wbn() && y().wbn() == o.y().wbn() && z().wbn() == o.z().wbn();
  }

  wide_jacobian_curve_point opposite() const {
    wide_jacobian_curve_point ret;
    ret.x() = x();
    ret.y() = gfp::mgry_opposite(y());
    ret.z() = z();
    return ret;
  }

  auto& x() { return x_; }
  auto& y() { return y_; }
  auto& z() { return z_; }

  auto const& x() const { return x_; }
  auto const& y() const { return y_; }
  auto const& z() const { return z_; }

private:
  WMBN x_;
  WMBN y_;
  WMBN z_;
};

} // ecsimd

#endif
