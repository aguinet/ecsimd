#ifndef FPSIMD_BIGINT_MUL_H
#define FPSIMD_BIGINT_MUL_H

#include <eve/wide.hpp>
#include <eve/function/all.hpp>
#include <eve/detail/meta.hpp>
#include <eve/traits/cardinal.hpp>

#include <ecsimd/bignum.h>

#include <immintrin.h>

#include <bit>
#include <cassert>

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

static auto sqr_wide(eve::wide<uint32_t, eve::fixed<4>> const a) {
#ifdef __AVX2__
  const __m128i a_ = std::bit_cast<__m128i>(a);
  const __m256i a64 = _mm256_cvtepu32_epi64(a);
  const __m256i ret = _mm256_mul_epu32(a64, a64);
  return std::bit_cast<eve::wide<uint64_t, eve::fixed<4>>>(ret);
#else
  return mul_wide(a,a);
#endif
}

static auto mullow(eve::wide<uint64_t, eve::fixed<4>> const a, eve::wide<uint64_t, eve::fixed<4>> const b) {
  const __m256i a_ = std::bit_cast<__m256i>(a);
  const __m256i b_ = std::bit_cast<__m256i>(b);
  const __m256i ret = _mm256_mul_epu32(a_, b_);
  return std::bit_cast<eve::wide<uint64_t, eve::fixed<4>>>(ret);
}

template <concepts::wide_bignum WBN>
static auto zext_u32x64(WBN const& v)
{
  using limb_type = bn_limb_t<WBN>;
  using cardinal = eve::cardinal_t<WBN>;
  constexpr size_t nlimbs = bn_nlimbs<WBN>;
  using half_limb_type = eve::detail::downgrade_t<limb_type>;

  using ret_type = eve::wide<bignum<limb_type, nlimbs*2>, cardinal>;

  const eve::wide<limb_type, cardinal> low_mask(std::numeric_limits<half_limb_type>::max());

  auto ret = eve::zero(eve::as<ret_type>());
  eve::detail::for_<0, 1, nlimbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    const auto l = eve::get<i>(v);
    eve::get<2*i>(ret)   = l & low_mask;
    eve::get<2*i+1>(ret) = l >> std::numeric_limits<half_limb_type>::digits;
  });
  return ret;
}

template <concepts::wide_bignum WBN>
static auto trunc_u64x32(WBN const& v)
{
  using limb_type = bn_limb_t<WBN>;
  using cardinal = eve::cardinal_t<WBN>;
  constexpr size_t nlimbs = bn_nlimbs<WBN>;
  static_assert(nlimbs % 2 == 0);
  using half_limb_type = eve::detail::downgrade_t<limb_type>;

  using wide_limb_type = eve::wide<limb_type, cardinal>;
  using ret_type = eve::wide<bignum<limb_type, nlimbs/2>, cardinal>;

  // TODO: generate (true,false)*(cardinal/2) at compile-time
  //const auto lm = eve::as_logical_t<>{true, false, true, false};

  auto ret = eve::zero(eve::as<ret_type>());
  eve::detail::for_<0, 1, nlimbs/2>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    const auto l0 = eve::get<2*i>(v);
    const auto l1 = eve::get<2*i+1>(v);

    // TODO: generalize
    const auto l0_ = std::bit_cast<__m256i>(l0);
    const auto l1_ = std::bit_cast<__m256i>(l1 << std::numeric_limits<half_limb_type>::digits);
    const auto l_  = _mm256_blend_epi32(l0_, l1_, 0b10101010);
    eve::get<i>(ret) = std::bit_cast<wide_limb_type>(l_);
  });
  return ret;
}

template <concepts::wide_bignum WBN>
static auto
  mul_u32_zext(WBN const& a, WBN const& b)
{
  using limb_type = bn_limb_t<WBN>;
  static_assert(std::is_same_v<limb_type, uint64_t>);
  using half_limb_type = eve::detail::downgrade_t<limb_type>;
  constexpr size_t nlimbs = bn_nlimbs<WBN>;
  constexpr auto half_limb_bits = std::numeric_limits<half_limb_type>::digits;
  using cardinal = eve::cardinal_t<WBN>;
  using ret_type = eve::wide<bignum<limb_type, nlimbs*2>, cardinal>;

  auto ret = eve::zero(eve::as<ret_type>());
  const eve::wide<limb_type, cardinal> low_mask(std::numeric_limits<half_limb_type>::max());

  eve::detail::for_<0, 1, nlimbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;
    auto highprev = eve::zero(eve::as<eve::wide<limb_type, cardinal>>());
    eve::detail::for_<0, 1, nlimbs>([&](auto j_) EVE_LAMBDA_FORCEINLINE {
      constexpr auto j = decltype(j_)::value;
      constexpr auto retlimb = i+j;

      auto t = mullow(eve::get<i>(a), eve::get<j>(b));
      t += eve::get<retlimb>(ret);
      t += highprev;
      eve::get<retlimb>(ret) = t & low_mask;

      highprev = t >> half_limb_bits;
    });
    eve::get<i + nlimbs>(ret) = highprev;
  });

  return ret;
}

template <concepts::wide_bignum WBN>
static auto
  mul(WBN const& a, WBN const& b)
{
  const auto a_half = zext_u32x64(a);
  const auto b_half = zext_u32x64(b);
  const auto m = mul_u32_zext(a_half,b_half);
  return trunc_u64x32(m);
}

template <concepts::wide_bignum WBN>
static auto
  square_u32_zext(WBN const& a)
{
  using limb_type = bn_limb_t<WBN>;
  using half_limb_type = eve::detail::downgrade_t<limb_type>;
  constexpr size_t nlimbs = bn_nlimbs<WBN>;
  constexpr auto half_limb_bits = std::numeric_limits<half_limb_type>::digits;
  using cardinal = eve::cardinal_t<WBN>;
  using WL = eve::wide<limb_type, cardinal>;
  using ret_type = eve::wide<bignum<limb_type, nlimbs*2>, cardinal>;

  auto ret = eve::zero(eve::as<ret_type>());
  const eve::wide<limb_type, cardinal> low_mask(std::numeric_limits<half_limb_type>::max());

  // Compute the cross products
  eve::detail::for_<0,1,nlimbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;

    const auto la = eve::get<i>(a);
    auto t = mullow(la, la);
    t += eve::get<2*i>(ret);
    eve::get<2*i>(ret) = t & low_mask;

    WL prevs[2];
    prevs[0] = t >> half_limb_bits;
    prevs[1] = eve::zero(eve::as<WL>());

    eve::detail::for_<i+1,1,nlimbs>([&](auto j_) EVE_LAMBDA_FORCEINLINE {
      constexpr auto j = decltype(j_)::value;
      static_assert(i != j);
      constexpr auto retlimb = i+j;

      auto t = mullow(eve::get<i>(a), eve::get<j>(a));
      const auto carry = t >> ((2*half_limb_bits)-1);
      t <<= 1;
      t += eve::get<retlimb>(ret);
      t += prevs[0];

      eve::get<retlimb>(ret) = t & low_mask;

      prevs[0] = prevs[1];
      prevs[0] += t >> half_limb_bits;
      prevs[1] = carry;
    });

    eve::get<i+nlimbs>(ret) += prevs[0]; // TODO: carry?
    if constexpr ((i+nlimbs+1) < (2*nlimbs)) {
      eve::get<i+nlimbs+1>(ret) = prevs[1];
    }
  });
  return ret;
}

template <concepts::wide_bignum WBN>
static auto
  square(WBN const& a)
{
  const auto a_half = zext_u32x64(a);
  const auto s = square_u32_zext(a_half);
  return trunc_u64x32(s);
}

template <concepts::wide_bignum WBN>
static auto
  limb_mul_zext(WBN const& a, eve::wide<bn_limb_t<WBN>, eve::cardinal_t<WBN>> const b)
{
  // b must contains half_limb_type constants zero extended to limb_type
  // TODO: assert for this
  using limb_type = bn_limb_t<WBN>;
  using half_limb_type = eve::detail::downgrade_t<limb_type>;
  constexpr size_t nlimbs = bn_nlimbs<WBN>;
  constexpr auto half_limb_bits = std::numeric_limits<half_limb_type>::digits;
  using cardinal = eve::cardinal_t<WBN>;
  using ret_type = eve::wide<bignum<limb_type, nlimbs+1>, cardinal>;

  auto ret = eve::zero(eve::as<ret_type>());
  const eve::wide<limb_type, cardinal> low_mask(std::numeric_limits<half_limb_type>::max());

  eve::wide<limb_type, cardinal> highprev;
  eve::detail::for_<0, 1, nlimbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
    constexpr auto i = decltype(i_)::value;

    auto t = mullow(eve::get<i>(a), b);
    if constexpr (i > 0) {
      t += highprev;
    }
    eve::get<i>(ret) = t & low_mask;
    highprev = t >> half_limb_bits;
  });
  eve::get<nlimbs>(ret) = highprev & low_mask;
  return ret;
}

template <concepts::wide_bignum WBN>
static auto
  limb_mul(WBN const& a, eve::wide<bn_limb_t<WBN>, eve::cardinal_t<WBN>> const b)
{
  using limb_type = bn_limb_t<WBN>;
  using dbl_limb_type = eve::detail::upgrade_t<limb_type>;

  // For testing purposes
  const auto azext = zext_u32x64(a);
  const auto rzext = limb_mul_zext(azext, b);
  return trunc_u64x32(pad<1>(rzext));
}

} // ecsmid

#endif
