#include <ecsimd/ops.h>
#include <ecsimd/serialization.h>
#include <ecsimd/literals.h>

#include <eve/function/all.hpp>

#include <gtest/gtest.h>

using namespace ecsimd;
using namespace ecsimd::literals;

template <class Bignum>
static auto wide_bignum_set1(std::array<uint8_t, sizeof(Bignum)> const& bytes) {
  const auto BN = bn_from_bytes_BE<Bignum>(bytes);
  return eve::wide<Bignum>{ [&](auto i, auto) { return BN; } };
}

template <class Bignum, size_t N>
static auto DoAdd(std::array<uint8_t, N> const& v0, std::array<uint8_t, N> const& v1)
{
  static_assert(N == sizeof(Bignum));
  const auto wv0 = wide_bignum_set1<bignum_128>(v0);
  const auto wv1 = wide_bignum_set1<bignum_128>(v1);

  return ecsimd::add(wv0, wv1);
}

template <class Bignum, size_t N>
static auto DoSub(std::array<uint8_t, N> const& v0, std::array<uint8_t, N> const& v1)
{
  static_assert(N == sizeof(Bignum));
  const auto wv0 = wide_bignum_set1<bignum_128>(v0);
  const auto wv1 = wide_bignum_set1<bignum_128>(v1);

  return ecsimd::sub(wv0, wv1);
}

TEST(Ops128, Add) {
  EXPECT_TRUE(eve::all(DoAdd<bignum_128>("00000000000000000000000500000005"_hex, "0000000000000000FFFFFFFFFFFFFFFF"_hex) ==
    wide_bignum_set1<bignum_128>("00000000000000010000000500000004"_hex)));
  EXPECT_TRUE(eve::all(DoAdd<bignum_128>("909680e1f399ca5916134a18b816399b"_hex, "0e36dfecf5e7f74363c453efc1cbc153"_hex) ==
    wide_bignum_set1<bignum_128>("9ecd60cee981c19c79d79e0879e1faee"_hex)));
}

TEST(Ops128, Sub) {
  EXPECT_TRUE(eve::all(DoSub<bignum_128>("00000000000000000000000500000005"_hex, "0000000000000000FFFFFFFFFFFFFFFF"_hex) ==
    wide_bignum_set1<bignum_128>("ffffffffffffffff0000000500000006"_hex)));
}
