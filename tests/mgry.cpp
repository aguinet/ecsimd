#include <ecsimd/mgry.h>
#include <ecsimd/mgry_mul.h>
#include <ecsimd/mgry_ops.h>
#include <ecsimd/serialization.h>
#include <ecsimd/literals.h>

#include <ctbignum/decimal_literals.hpp>
#include <ctbignum/montgomery.hpp>
#include <ctbignum/io.hpp>
#include <ctbignum/relational_ops.hpp>

#include <eve/function/all.hpp>

#include <gtest/gtest.h>

#include "tests.h"

using namespace ecsimd;
using namespace ecsimd::literals;

using namespace cbn::literals;

namespace {
struct P {
  static constexpr auto value = bn_from_bytes_BE<bignum_256>("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F"_hex);
};
constexpr auto P_cbn = P::value.cbn();
using P_cbn_is = array_to_integer_sequence_t<P_cbn>;
}

TEST(Mgry, FromTo) {
  using WBN = wide_bignum<bignum_256>;
  using WMBN = wide_mgry_bignum<WBN, P>;
  constexpr std::array<uint8_t, 32> bns[] = {
    "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"_hex,
    "0168db3a8eca3fd7d4d08943182e189aef318068ba8853d77cb49c17bae00c0e"_hex,
    "2714dac0b974321b75d6ef64e7c3b118adb2801bf674282df5712cd2af390f79"_hex,
    "a3fc64fece6f3e1effab4045a9a54faa49a228f787025f0ecb761145755cb2d0"_hex,
    "3af178b78710adae9cc096188ed09c210078aaa7e965ef83d22a91f21fec4eb5"_hex,
    "688c743cde3987e299d2b028038ddc12dc02e7033c9d3c8f4d20edf9544232aa"_hex,
    "45e29166c6441f0fd27e3b85a205f1e102b025cc8e8ea158ab4885a22ed68905"_hex,
  };
  for (auto const& bn: bns) {
    const auto a = wide_bignum_set1<WBN>(bn);
    const auto ma = WMBN::from_classical(a);
    const auto maback = ma.to_classical();
    EXPECT_TRUE(eve::all(maback == a));
  }
}

static void TestMul(bignum_256 const& a, bignum_256 const& b)
{
  using WBN = wide_bignum<bignum_256>;
  const auto vvec = ecsimd::details::mgry_mul<P>(WBN{a}, WBN{b});

  const auto a_cbn = a.cbn();
  const auto b_cbn = b.cbn();
  const auto v = cbn::montgomery_mul(a_cbn,b_cbn,P_cbn_is{});

  EXPECT_EQ(vvec.get(0).cbn(), v);
}

TEST(Mgry, Mul) {
  {
    const auto a = bn_from_bytes_BE<bignum_256>("0000000000000000000000000000000000000000000000000000000000000004"_hex);
    const auto b = bn_from_bytes_BE<bignum_256>("0000000000000000000000000000000000000000000000000000000000000005"_hex);
    TestMul(a, b);
  }

  {
    const auto a = bn_from_bytes_BE<bignum_256>("00000000000AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex);
    const auto b = bn_from_bytes_BE<bignum_256>("00000000000BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"_hex);
    TestMul(a, b);
  }

  {
    const auto a = bn_from_bytes_BE<bignum_256>("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2E"_hex);
    const auto b = bn_from_bytes_BE<bignum_256>("0000000000000000000000000000000000000000000000000000000000000001"_hex);
    TestMul(a, b);
  }

  {
    const auto a = bn_from_bytes_BE<bignum_256>("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2E"_hex);
    const auto b = bn_from_bytes_BE<bignum_256>("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2D"_hex);
    TestMul(a, b);
  }

  {
    const auto a = bn_from_bytes_BE<bignum_256>("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2E"_hex);
    const auto b = bn_from_bytes_BE<bignum_256>("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2E"_hex);
    TestMul(a, b);
  }
}

TEST(Mgry, Ops) {
  using WBN = wide_bignum<bignum_256>;
  using WMBN = wide_mgry_bignum<WBN, P>;

  const auto a = wide_bignum_set1<WBN>("FFFFFFFFFFFFFFFFFFFFFF000000000000000000000000000000000000000004"_hex);
  const auto b = wide_bignum_set1<WBN>("FFFFFFFFFFFFFFFFFFFFFF000000000000000000000000000000000000000005"_hex);

  const auto ma = WMBN::from_classical(a);
  const auto mb = WMBN::from_classical(b);

  const auto add = mgry_add(ma, mb);
  EXPECT_TRUE(eve::all(add.to_classical() == wide_bignum_set1<WBN>("fffffffffffffffffffffe0000000000000000000000000000000001000003da"_hex)));

  auto sub = mgry_sub(ma, mb);
  EXPECT_TRUE(eve::all(sub.to_classical() == wide_bignum_set1<WBN>("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2e"_hex)));

  sub = mgry_sub(mb, ma);
  EXPECT_TRUE(eve::all(sub.to_classical() == wide_bignum_set1<WBN>("0000000000000000000000000000000000000000000000000000000000000001"_hex)));
}
