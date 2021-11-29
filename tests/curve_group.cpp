#include <ecsimd/curve_point.h>
#include <ecsimd/jacobian_curve_point.h>
#include <ecsimd/curve_group.h>
#include <ecsimd/curve_nist_p256.h>
#include <ecsimd/serialization.h>
#include <ecsimd/swap.h>
#include <ecsimd/literals.h>

#include <eve/function/all.hpp>

#include <gtest/gtest.h>

#include "tests.h"

using namespace ecsimd;
using namespace ecsimd::literals;

template <class WJCP>
static void dump_pt(WJCP const& pt) {
  const auto x = bn_to_bytes_BE(pt.x().to_classical().get(0));
  const auto y = bn_to_bytes_BE(pt.y().to_classical().get(0));
  const auto z = bn_to_bytes_BE(pt.z().to_classical().get(0));

  for (uint8_t v: x) {
    printf("%02X", v);
  }
  printf(" ");
  for (uint8_t v: y) {
    printf("%02X", v);
  }
  printf(" ");
  for (uint8_t v: z) {
    printf("%02X", v);
  }
  puts("");
}

TEST(CurveGroup, DBLU) {
  using Curve = curve_nist_p256;
  using CurveGroup = curve_group<Curve>;
  using WBN = curve_wide_bn_t<Curve>;

  auto WJG = CurveGroup::WJG();
  const auto WJdblG = CurveGroup::DBLU(WJG);
  EXPECT_TRUE(eve::all(WJG.z().wbn() == WJdblG.z().wbn()));
  EXPECT_TRUE(eve::all(WJG.to_affine() == CurveGroup::WG()));

  const auto WdblG = WJdblG.to_affine();
  // 2*G in affine coordinates
  EXPECT_TRUE(eve::all(WdblG.x() == wide_bignum_set1<WBN>("7cf27b188d034f7e8a52380304b51ac3c08969e277f21b35a60b48fc47669978"_hex)));
  EXPECT_TRUE(eve::all(WdblG.y() == wide_bignum_set1<WBN>("07775510db8ed040293d9ac69f7430dbba7dade63ce982299e04b79d227873d1"_hex)));
}

TEST(CurveGroup, ZADDU) {
  using Curve = curve_nist_p256;
  using CurveGroup = curve_group<Curve>;
  using WBN = curve_wide_bn_t<Curve>;

  auto WJG = CurveGroup::WJG();
  const auto WJdblG = CurveGroup::DBLU(WJG);
  const auto WJP = CurveGroup::ZADDU(WJG, WJdblG);
  EXPECT_TRUE(eve::all(WJP.z().wbn() == WJG.z().wbn()));

  auto WP = WJP.to_affine();
  // 3*G in affine coordinates
  const auto W3G_x = wide_bignum_set1<WBN>("5ecbe4d1a6330a44c8f7ef951d4bf165e6c6b721efada985fb41661bc6e7fd6c"_hex);
  const auto W3G_y = wide_bignum_set1<WBN>("8734640c4998ff7e374b06ce1a64a2ecd82ab036384fb83d9a79b127a27d5032"_hex);
  EXPECT_TRUE(eve::all(WP.x() == W3G_x));
  EXPECT_TRUE(eve::all(WP.y() == W3G_y));

  WJG = CurveGroup::WJG();
  const auto WJ3P = CurveGroup::TRPLU(WJG);
  WP = WJ3P.to_affine();
  EXPECT_TRUE(eve::all(WP.x() == W3G_x));
  EXPECT_TRUE(eve::all(WP.y() == W3G_y));
}

TEST(CurveGroup, ZDAU) {
  using Curve = curve_nist_p256;
  using CurveGroup = curve_group<Curve>;
  using WBN = curve_wide_bn_t<Curve>;

  auto WJG = CurveGroup::WJG();
  auto WdblG = CurveGroup::DBLU(WJG);
  const auto WJP = CurveGroup::ZDAU(WdblG, WJG);
  EXPECT_TRUE(eve::all(WJP.z().wbn() == WJG.z().wbn()));

  auto WP = WJP.to_affine();
  // 5*G in affine coordinates
  const auto W5G_x = wide_bignum_set1<WBN>("51590b7a515140d2d784c85608668fdfef8c82fd1f5be52421554a0dc3d033ed"_hex);
  const auto W5G_y = wide_bignum_set1<WBN>("e0c17da8904a727d8ae1bf36bf8a79260d012f00d4d80888d1d0bb44fda16da4"_hex);
  EXPECT_TRUE(eve::all(WP.x() == W5G_x));
  EXPECT_TRUE(eve::all(WP.y() == W5G_y));
}

TEST(CurveGroup, Swap) {
  using Curve = curve_nist_p256;
  using CurveGroup = curve_group<Curve>;
  using WBN = curve_wide_bn_t<Curve>;

  auto WJG = CurveGroup::WJG();
  auto WdblG = CurveGroup::DBLU(WJG);

  const auto Z = eve::zero(eve::as<cmp_res_t<WBN>>());

  auto a = WJG;
  auto b = WdblG;
  swap_if(Z, a, b);
  EXPECT_TRUE(eve::all(a == WJG));
  EXPECT_TRUE(eve::all(b == WdblG));

  swap_if(!Z, a, b);
  EXPECT_TRUE(eve::all(a == WdblG));
  EXPECT_TRUE(eve::all(b == WJG));
}

TEST(CurveGroup, ScalarMult) {
  using Curve = curve_nist_p256;
  using CurveGroup = curve_group<Curve>;
  using WBN = curve_wide_bn_t<Curve>;
  using BN = typename WBN::value_type;

  auto WJG = CurveGroup::WJG();
  {
    const auto xs = bn_from_bytes_BE<BN>("0000000000000000000000000000000000000000000000000000000000000005"_hex);
    const auto x = WBN{xs};
    const auto WJP = CurveGroup::scalar_mult(x, WJG);
    const auto WP = WJP.to_affine();

    const auto WJPs = CurveGroup::scalar_mult_1s(xs, WJG);
    const auto WPs = WJPs.to_affine();

    // 5*G in affine coordinates
    const auto W5G_x = wide_bignum_set1<WBN>("51590b7a515140d2d784c85608668fdfef8c82fd1f5be52421554a0dc3d033ed"_hex);
    const auto W5G_y = wide_bignum_set1<WBN>("e0c17da8904a727d8ae1bf36bf8a79260d012f00d4d80888d1d0bb44fda16da4"_hex);
    EXPECT_TRUE(eve::all(WP.x() == W5G_x));
    EXPECT_TRUE(eve::all(WP.y() == W5G_y));
    EXPECT_TRUE(eve::all(WPs.x() == W5G_x));
    EXPECT_TRUE(eve::all(WPs.y() == W5G_y));
  }
  {
    const auto xs = bn_from_bytes_BE<BN>("0bc1b1f28709decb543d9677d2cc9942348f6b984deff409430740942ff38827"_hex);
    const auto x = WBN{xs};
    const auto WJP = CurveGroup::scalar_mult(x, WJG);
    const auto WP = WJP.to_affine();

    const auto WJPs = CurveGroup::scalar_mult_1s(xs, WJG);
    const auto WPs = WJPs.to_affine();

    const auto WP_x = wide_bignum_set1<WBN>("1b7721565b2c4a9f203bbccc6b531df2789fde0d135c76db71e4a7bbab9e85b2"_hex);
    const auto WP_y = wide_bignum_set1<WBN>("393655bcc30f67f3a4e257b39685657d7c8df7b2a132b49c848003e300c8dcd1"_hex);
    EXPECT_TRUE(eve::all(WP.x() == WP_x));
    EXPECT_TRUE(eve::all(WP.y() == WP_y));
    EXPECT_TRUE(eve::all(WPs.x() == WP_x));
    EXPECT_TRUE(eve::all(WPs.y() == WP_y));
  }
  {
    const auto xs = bn_from_bytes_BE<BN>("0a891cecc2bf13b0aca744434a9c9f4bd7bf5c8ed86e2f76e7df72bad813bd80"_hex);
    const auto x = WBN{xs};
    const auto WJP = CurveGroup::scalar_mult(x, WJG);
    const auto WP = WJP.to_affine();

    const auto WJPs = CurveGroup::scalar_mult_1s(xs, WJG);
    const auto WPs = WJPs.to_affine();

    const auto WP_x = wide_bignum_set1<WBN>("f411d79e2997b2954975046d23b0e4a69ce580a4a81e1bed18fef6fd9ea4a912"_hex);
    const auto WP_y = wide_bignum_set1<WBN>("43895f527937e816c3d7c0a2370002796d3cd4860cb034df86cbe7da227d9113"_hex);
    EXPECT_TRUE(eve::all(WP.x() == WP_x));
    EXPECT_TRUE(eve::all(WP.y() == WP_y));
    EXPECT_TRUE(eve::all(WPs.x() == WP_x));
    EXPECT_TRUE(eve::all(WPs.y() == WP_y));
  }
}
