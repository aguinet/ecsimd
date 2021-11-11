#include <ecsimd/add.h>
#include <ecsimd/sub.h>
#include <ecsimd/mul.h>
#include <ecsimd/cmp.h>
#include <ecsimd/modular.h>
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
static auto DoMul = [](auto const& v0, auto const& v1) { return DoFunc<WBN>(v0, v1, &ecsimd::mul<bn_nlimbs<WBN>>); };
template <concepts::wide_bignum WBN>
static auto DoSub = [](auto const& v0, auto const& v1) { return DoFunc<WBN>(v0, v1, &ecsimd::sub_no_carry<WBN>); };
template <concepts::wide_bignum WBN>
static auto DoSubIfAbove = [](auto const& v0, auto const& v1) { return DoFunc<WBN>(v0, v1, &ecsimd::sub_if_above<WBN>); };
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
  using WDBN = wide_bignum<bignum_256>;

  // Additions
  EXPECT_TRUE(eve::all(DoAdd<WBN>("00000000000000000000000500000005"_hex, "0000000000000000FFFFFFFFFFFFFFFF"_hex) ==
    wide_bignum_set1<WBN>("00000000000000010000000500000004"_hex)));
  EXPECT_TRUE(eve::all(DoAdd<WBN>("909680e1f399ca5916134a18b816399b"_hex, "0e36dfecf5e7f74363c453efc1cbc153"_hex) ==
    wide_bignum_set1<WBN>("9ecd60cee981c19c79d79e0879e1faee"_hex)));

  // Substractions
  EXPECT_TRUE(eve::all(DoSub<WBN>("00000000000000000000000500000005"_hex, "0000000000000000FFFFFFFFFFFFFFFF"_hex) ==
    wide_bignum_set1<WBN>("ffffffffffffffff0000000500000006"_hex)));

  // Substracte if above
  EXPECT_TRUE(eve::all(DoSubIfAbove<WBN>("F0000000000000000000000000000005"_hex, "F0000000000000000000000000000004"_hex) ==
    wide_bignum_set1<WBN>("00000000000000000000000000000001"_hex)));
  EXPECT_TRUE(eve::all(DoSubIfAbove<WBN>("F0000000000000000000000000000004"_hex, "F0000000000000000000000000000004"_hex) ==
    wide_bignum_set1<WBN>("00000000000000000000000000000000"_hex)));
  EXPECT_TRUE(eve::all(DoSubIfAbove<WBN>("F0000000000000000000000000000003"_hex, "F0000000000000000000000000000004"_hex) ==
    wide_bignum_set1<WBN>("F0000000000000000000000000000003"_hex)));

  {
    constexpr std::array<uint8_t, 16> a[] = {
      "F0000000000000000000000000000005"_hex,
      "F0000000000000000000000000000004"_hex,
      "F0000000000000000000000000000003"_hex,
      "F0000000000000000000000000000002"_hex
    };
    const WBN wa{[&](auto i, auto _) { return bn_from_bytes_BE<bignum_128>(a[i]); }};

    constexpr std::array<uint8_t, 16> b[] = {
      "F0000000000000000000000000000004"_hex,
      "F0000000000000000000000000000004"_hex,
      "F0000000000000000000000000000004"_hex,
      "F0000000000000000000000000000004"_hex
    };
    const WBN wb{[&](auto i, auto _) { return bn_from_bytes_BE<bignum_128>(b[i]); }};

    constexpr std::array<uint8_t, 16> res[] = {
      "00000000000000000000000000000001"_hex,
      "00000000000000000000000000000000"_hex,
      "F0000000000000000000000000000003"_hex,
      "F0000000000000000000000000000002"_hex
    };
    const WBN wres{[&](auto i, auto _) { return bn_from_bytes_BE<bignum_128>(res[i]); }};

    EXPECT_TRUE(eve::all(sub_if_above(wa, wb) == wres));
  }

  // Multiplications
  EXPECT_TRUE(eve::all(DoMul<WBN>("ffffffffffffffffffffffffffffffff"_hex, "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"_hex) ==
    wide_bignum_set1<WDBN>("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEED11111111111111111111111111111112"_hex)));

  // Comparaisons
  EXPECT_TRUE(eve::all(DoLT<WBN> ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex, "BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex)));
  EXPECT_TRUE(eve::all(DoLTE<WBN>("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex, "BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex)));
  EXPECT_TRUE(eve::all(DoLTE<WBN>("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex)));
  EXPECT_TRUE(eve::all(DoGT<WBN> ("BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex)));
  EXPECT_TRUE(eve::all(DoGTE<WBN>("BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex)));
  EXPECT_TRUE(eve::all(DoGTE<WBN>("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex)));

  EXPECT_TRUE(eve::none(DoLT<WBN> ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB"_hex, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_hex)));
}

TEST(Ops256, Binops) {
  using WBN = wide_bignum<bignum_256>;
  using WDBN = wide_bignum<bignum_512>;

  // Multiplications
  EXPECT_TRUE(eve::all(DoMul<WBN>(
      "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"_hex,
      "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"_hex) ==
    wide_bignum_set1<WDBN>("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEED1111111111111111111111111111111111111111111111111111111111111112"_hex)));
}

TEST(Ops256, Mod) {
  using WBN = wide_bignum<bignum_256>;
  const auto p = wide_bignum_set1<WBN>("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F"_hex);

  {
    const auto a = wide_bignum_set1<WBN>("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2E"_hex);
    const auto b = wide_bignum_set1<WBN>("0000000000000000000000000000000000000000000000000000000000000002"_hex);
    EXPECT_TRUE(eve::all(mod_add(a,b,p) == wide_bignum_set1<WBN>("0000000000000000000000000000000000000000000000000000000000000001"_hex)));
  }

  {
    const auto a = wide_bignum_set1<WBN>("fffffffffffffffffffffffffffffffffffffffffffffffffffffff000000000"_hex);
    const auto b = wide_bignum_set1<WBN>("ffeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"_hex);
    const auto sum = mod_add(a,b,p);

    {
      const auto buf = bn_to_bytes_BE(sum.get(0));
      for (uint8_t v: buf) {
        printf("%02X", v);
      }
      printf("\n");
    }
    EXPECT_TRUE(eve::all(sum == wide_bignum_set1<WBN>("ffeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeedfeeeef2bf"_hex)));
  }
}
