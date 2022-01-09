#include <ecsimd/curve_group.h>
#include <ecsimd/curve_nist_p256.h>

using namespace ecsimd;

using Curve = curve_nist_p256;
using CurveGroup = curve_group<Curve>;
using WBN = curve_wide_bn_t<Curve>;
using WJCP = wide_jacobian_curve_point<Curve>;

auto scalar_mult_p256(WBN const& x, WJCP const& P) {
  return curve_group<Curve>::scalar_mult(x, P);
}
