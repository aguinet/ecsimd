#include <ecsimd/mgry.h>
#include <ecsimd/mgry_mul.h>
#include <ecsimd/serialization.h>
#include <ecsimd/literals.h>

#include <ctbignum/decimal_literals.hpp>
#include <ctbignum/montgomery.hpp>
#include <ctbignum/io.hpp>
#include <ctbignum/relational_ops.hpp>

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
/*
TEST(Mgry, FromTo) {
  const auto a = wide_bignum_set1<bignum_256, eve::fixed<4>>("eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"_hex);
  const auto ma = wide_mgry_bignum<wide_bignum<bignum_256>, P>::from_classical(a);
}
*/

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
