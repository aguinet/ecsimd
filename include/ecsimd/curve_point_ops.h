#ifndef ECSIMD_CURVE_POINT_OPS_H
#define ECSIMD_CURVE_POINT_OPS_H

#include <ecsimd/curve.h>
#include <ecsimd/curve_point.h>
#include <ecsimd/curve_group.h>

#include <optional>

namespace ecsimd {

template <concepts::curve Curve>
std::optional<wide_curve_point<Curve>> wide_curve_point<Curve>::from_x(typename wide_curve_point<Curve>::WBN const& x) {
  const auto y = curve_group<Curve>::compute_y(x);
  if (!y) {
    return {};
  }
  wide_curve_point ret;
  ret.x() = x;
  ret.y() = *y;
  return {ret};
}

} // ecsimd

#endif
