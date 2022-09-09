#include <ecsimd/curve_point.h>
#include <ecsimd/curve_point_ops.h>
#include <ecsimd/jacobian_curve_point.h>
#include <ecsimd/curve_nist_p256.h>
#include <ecsimd/serialization.h>
#include <ecsimd/literals.h>

#include <eve/module/core/regular/all.hpp>

#include <gtest/gtest.h>

#include "tests.h"

using namespace ecsimd;
using namespace ecsimd::literals;

TEST(CurvePoint, FromX) {
  using Curve = curve_nist_p256;
  using WideCurvePoint = wide_curve_point<curve_nist_p256>;
  using WBN = curve_wide_bn_t<Curve>;

  const auto x = wide_bignum_set1<WBN>("ce11d601ec0e947529e66021a0cd3d57518d58d0d5f2eb7ed75805d78c986e60"_hex);
  const auto opt = WideCurvePoint::from_x(x);
  EXPECT_TRUE(opt.has_value());
  EXPECT_TRUE(eve::all(opt->y() == wide_bignum_set1<WBN>("f2a40cfbb248ae2c7749c76641b51b7137ccad8916931adf83b857e418fad591"_hex)));
}

TEST(JacobianCurvePoint, ToFromAffine) {
  using Curve = curve_nist_p256;
  using WideCurvePoint = wide_curve_point<curve_nist_p256>;
  using WideJacobianCurvePoint = wide_jacobian_curve_point<curve_nist_p256>;
  using WBN = curve_wide_bn_t<Curve>;

  const auto x = wide_bignum_set1<WBN>("ce11d601ec0e947529e66021a0cd3d57518d58d0d5f2eb7ed75805d78c986e60"_hex);
  const auto opt = WideCurvePoint::from_x(x);
  EXPECT_TRUE(opt.has_value());

  const auto jpt = WideJacobianCurvePoint::from_affine(*opt);
  const auto apt = jpt.to_affine();

  EXPECT_TRUE(eve::all(*opt == apt));
}
