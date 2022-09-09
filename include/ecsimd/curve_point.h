#ifndef ECSIMD_CURVE_POINT_H
#define ECSIMD_CURVE_POINT_H

#include <ecsimd/bignum.h>
#include <ecsimd/curve.h>

#include <optional>

namespace ecsimd {

template <concepts::curve Curve>
struct wide_curve_point
{
  using curve_type = Curve;
  using bignum_type = typename curve_type::bn_type;
  using WBN = wide_bignum<bignum_type>;

  wide_curve_point() = default;
  wide_curve_point(wide_curve_point const&) = default;

  wide_curve_point(WBN const& x, WBN const& y):
    x_(x),
    y_(y)
  { }

  static std::optional<wide_curve_point> from_x(WBN const& x);

  WBN const& x() const { return x_; }
  WBN const& y() const { return y_; }

  WBN& x() { return x_; }
  WBN& y() { return y_; }

  auto operator==(wide_curve_point const& o) const {
    return x() == o.x() && y() == o.y();
  }

private:
  WBN x_;
  WBN y_;
};

} // ecsimd

#endif
