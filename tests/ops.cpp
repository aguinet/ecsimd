#include <ecsimd/add.h>
#include <ecsimd/sub.h>
#include <ecsimd/mul.h>
#include <ecsimd/cmp.h>
#include <ecsimd/serialization.h>
#include <ecsimd/literals.h>

#include <eve/function/all.hpp>
#include <eve/function/none.hpp>

#include <gtest/gtest.h>

#include "tests.h"

using namespace ecsimd;
using namespace ecsimd::literals;

template <concepts::wide_bignum WBN, size_t N, class Func>
static auto DoFunc(std::array<uint8_t, N> const& v0, std::array<uint8_t, N> const& v1, Func F)
{
  static_assert(N == sizeof(typename WBN::value_type));
  const auto wv0 = wide_bignum_set1<WBN>(v0);
  const auto wv1 = wide_bignum_set1<WBN>(v1);

  return F(wv0, wv1);
}

template <concepts::wide_bignum WBN>
static auto DoAdd = [](auto const& v0, auto const& v1) { return DoFunc<WBN>(v0, v1, &ecsimd::add_no_carry<WBN>); };
template <concepts::wide_bignum WBN>
static auto DoSub = [](auto const& v0, auto const& v1) { return DoFunc<WBN>(v0, v1, &ecsimd::sub_no_carry<WBN>); };
template <concepts::wide_bignum WBN>
static auto DoLT  = [](auto const& v0, auto const& v1) { return DoFunc<WBN>(v0, v1, [](auto a, auto b) { return a < b; }); };
template <concepts::wide_bignum WBN>
static auto DoLTE = [](auto const& v0, auto const& v1) { return DoFunc<WBN>(v0, v1, [](auto a, auto b) { return a <= b; }); };
template <concepts::wide_bignum WBN>
static auto DoGT  = [](auto const& v0, auto const& v1) { return DoFunc<WBN>(v0, v1, [](auto a, auto b) { return a > b; }); };
template <concepts::wide_bignum WBN>
static auto DoGTE = [](auto const& v0, auto const& v1) { return DoFunc<WBN>(v0, v1, [](auto a, auto b) { return a >= b; }); };

TEST(Ops128, Binops) {
  using WBN = wide_bignum<bignum_128>;

  // Additions
  EXPECT_TRUE(eve::all(DoAdd<WBN>("00000000000000000000000500000005"_hex, "0000000000000000FFFFFFFFFFFFFFFF"_hex) ==
    wide_bignum_set1<WBN>("00000000000000010000000500000004"_hex)));
  EXPECT_TRUE(eve::all(DoAdd<WBN>("909680e1f399ca5916134a18b816399b"_hex, "0e36dfecf5e7f74363c453efc1cbc153"_hex) ==
    wide_bignum_set1<WBN>("9ecd60cee981c19c79d79e0879e1faee"_hex)));

  // Substractions
  EXPECT_TRUE(eve::all(DoSub<WBN>("00000000000000000000000500000005"_hex, "0000000000000000FFFFFFFFFFFFFFFF"_hex) ==
    wide_bignum_set1<WBN>("ffffffffffffffff0000000500000006"_hex)));

  // Comparaisons
  EXPECT_TRUE(eve::all(DoLT<WBN> ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex, "BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex)));
  EXPECT_TRUE(eve::all(DoLTE<WBN>("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex, "BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex)));
  EXPECT_TRUE(eve::all(DoLTE<WBN>("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex)));
  EXPECT_TRUE(eve::all(DoGT<WBN> ("BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex)));
  EXPECT_TRUE(eve::all(DoGTE<WBN>("BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex)));
  EXPECT_TRUE(eve::all(DoGTE<WBN>("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex)));

  EXPECT_TRUE(eve::none(DoLT<WBN> ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB"_hex, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex)));
}

TEST(Ops128, Mul) {
  using WBN = wide_bignum<bignum_128>;
  const auto a = wide_bignum_set1<WBN>("ffffffffffffffffffffffffffffffff"_hex);
  const auto b = wide_bignum_set1<WBN>("eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"_hex);

  {
    const auto res = ecsimd::mul(a,b);
    const auto buf = bn_to_bytes_BE(res.get(0));
    for (uint8_t v: buf) {
      printf("%02X", v);
    }
    printf("\n");
  }
}

TEST(Ops256, Mul) {
  using WBN = wide_bignum<bignum_256>;
  const auto a = wide_bignum_set1<WBN>("ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"_hex);
  const auto b = wide_bignum_set1<WBN>("eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"_hex);

  {
    const auto res = ecsimd::mul(a,b);
    const auto buf = bn_to_bytes_BE(res.get(0));
    for (uint8_t v: buf) {
      printf("%02X", v);
    }
    printf("\n");
  }
}
