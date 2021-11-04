#ifndef FPSIMD_BIGINT_MUL_H
#define FPSIMD_BIGINT_MUL_H

#include <eve/wide.hpp>
#include <ecsimd/bignum.h>

#include <immintrin.h>

#include <bit>

namespace ecsimd {

static auto mul_wide_raw_m128(eve::wide<uint32_t, eve::fixed<4>> const a, eve::wide<uint32_t, eve::fixed<4>> const b) {
	const __m128i a_ = std::bit_cast<__m128i>(a);
	const __m128i b_ = std::bit_cast<__m128i>(b);
	const auto mul0 = _mm_mul_epu32(a, b);
	const auto mul1 = _mm_mul_epu32(
			_mm_shuffle_epi32(a,0b11110101),
			_mm_shuffle_epi32(b,0b11110101));
  return std::make_tuple(mul0, mul1);
}

static auto mul_wide(eve::wide<uint32_t, eve::fixed<4>> const a, eve::wide<uint32_t, eve::fixed<4>> const b) {
#ifdef __AVX2__
  const __m128i a_ = std::bit_cast<__m128i>(a);
  const __m128i b_ = std::bit_cast<__m128i>(b);
  const __m256i a64 = _mm256_cvtepu32_epi64(a);
  const __m256i b64 = _mm256_cvtepu32_epi64(b);
  const __m256i ret = _mm256_mul_epu32(a64, b64);
#else
  const auto [mul0, mul1] = mul_wide_raw_m128(a, b);
  struct { __m128i v0; __m128i v1;} ret = {
    _mm_unpacklo_epi64(mul0, mul1),
    _mm_unpackhi_epi64(mul0, mul1)
  };
#endif
  return std::bit_cast<eve::wide<uint64_t, eve::fixed<4>>>(ret);
}

// HACK: temporarely home-made for Nx32 bits integers.
template <size_t N>
__attribute__((noinline)) static auto
  mul(eve::wide<bignum<uint32_t, N>, eve::fixed<4>> const& a, eve::wide<bignum<uint32_t, N>, eve::fixed<4>> const& b)
{
  using limb_type = uint32_t;
  using dbl_limb_type = uint64_t;
  constexpr size_t nlimbs = N;
  constexpr auto limb_bits = std::numeric_limits<limb_type>::digits;
  using cardinal = eve::fixed<4>;
  using ret_type = eve::wide<bignum<uint32_t, N*2>, cardinal>;

  auto ret = eve::zero(eve::as<ret_type>());

  eve::detail::for_<0, 1, nlimbs>([&](auto i_) {
    constexpr auto i = decltype(i_)::value;
    auto highprev = eve::zero(eve::as<eve::wide<dbl_limb_type, cardinal>>());
    eve::detail::for_<0, 1, nlimbs>([&](auto j_) {
      constexpr auto j = decltype(j_)::value;
      constexpr auto retlimb = i+j;

      auto t = mul_wide(eve::get<i>(a), eve::get<j>(b));
      t += eve::convert(eve::get<retlimb>(ret), eve::as<dbl_limb_type>());
      t += highprev;
      // Chances are we can "hack" this to make it faster using an AVX2 shuffle + addition.
      eve::get<retlimb>(ret) = eve::convert(t, eve::as<limb_type>());

      highprev = t >> limb_bits;
    });
    eve::get<i + nlimbs>(ret) = eve::convert(highprev, eve::as<limb_type>());
  });
  return ret;
}

} // ecsmid

#endif
