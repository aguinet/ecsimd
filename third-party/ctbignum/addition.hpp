//
// This file is part of
//
// CTBignum 	
//
// C++ Library for Compile-Time and Run-Time Multi-Precision and Modular Arithmetic
// 
//
// This file is distributed under the Apache License, Version 2.0. See the LICENSE
// file for details.

#ifndef CT_ADDITION_HPP
#define CT_ADDITION_HPP

#include <ctbignum/bigint.hpp>
#include <ctbignum/config.hpp>
#include <ctbignum/slicing.hpp>

#include <algorithm>
#include <cstddef>
#include <bit>

namespace cbn {

template <typename T, size_t M, size_t N>
CBN_ALWAYS_INLINE
constexpr auto add(big_int<M, T> a, big_int<N, T> b) {
  constexpr auto L = std::max(M, N);
  return add_same(detail::pad<L - M>(a), detail::pad<L - N>(b));
}

template <typename T, size_t M, size_t N>
CBN_ALWAYS_INLINE
constexpr auto subtract(big_int<M, T> a, big_int<N, T> b) {
  constexpr auto L = std::max(M, N);
  return subtract_same(detail::pad<L - M>(a), detail::pad<L - N>(b));
}


template <typename T, size_t N>
CBN_ALWAYS_INLINE
constexpr auto add_same(big_int<N, T> a, big_int<N, T> b) {
  T carry{};
  big_int<N + 1, T> r{};

  for (auto i = 0U; i < N; ++i) {
    auto aa = a[i];
    auto sum = aa + b[i];
    auto res = sum + carry;
    carry = (sum < aa) | (res < sum);
    r[i] = res;
  }

  r[N] = carry;
  return r;
}

template <typename T, size_t N>
CBN_ALWAYS_INLINE
constexpr auto subtract_same(big_int<N, T> a, big_int<N, T> b) {
  T carry{};
  big_int<N + 1, T> r{};

  for (auto i = 0U; i < N; ++i) {
    auto aa = a[i];
    auto diff = aa - b[i];
    auto res = diff - carry;
    carry = (diff > aa) | (res > diff);
    r[i] = res;
  }

  r[N] = carry * static_cast<T>(-1); // sign extension
  return r;
}

template <typename T, size_t N>
CBN_ALWAYS_INLINE
constexpr auto add_ignore_carry(big_int<N, T> a, big_int<N, T> b) {
  T carry{};
  big_int<N, T> r{};

  for (auto i = 0U; i < N; ++i) {
    T aa = a[i];
    T sum = aa + b[i];
    T res = sum + carry;
    carry = (sum < aa) | (res < sum);
    r[i] = res;
  }

  return r;
}

template <typename T, size_t N>
constexpr auto subtract_ignore_carry(big_int<N, T> a, big_int<N, T> b) {
  T carry{};
  big_int<N, T> r{};

  for (auto i = 0U; i < N; ++i) {
    auto aa = a[i];
    auto diff = aa - b[i];
    auto res = diff - carry;
    carry = (diff > aa) | (res > diff);
    r[i] = res;
  }

  return r;
}

template <typename T, size_t N>
constexpr auto mod_add(big_int<N, T> a, big_int<N, T> b,
                       big_int<N, T> modulus) {
#ifdef __clang__
  if constexpr (std::is_same_v<T, uint64_t> && N == 4) {
    using u256 = _ExtInt(256);
    using u256_ar = big_int<N,T>;
    uint64_t carryin=0, carryout;
    uint64_t const* x = &a[0];
    uint64_t const* y = &b[0];
    u256_ar z;
    z[0] = __builtin_addcl(x[0], y[0], carryin, &carryout);
    carryin = carryout;
    z[1] = __builtin_addcl(x[1], y[1], carryin, &carryout);
    carryin = carryout;
    z[2] = __builtin_addcl(x[2], y[2], carryin, &carryout);
    carryin = carryout;
    z[3] = __builtin_addcl(x[3], y[3], carryin, &carryout);

    const uint64_t mask = -carryout;
    const u256_ar m256{mask, mask, mask, mask};
    const u256 pm = std::bit_cast<u256>(modulus) & std::bit_cast<u256>(m256);
    const auto retmp = std::bit_cast<u256_ar>(std::bit_cast<u256>(z)-pm);
    return retmp;
  }
#endif
  T carry{};
  big_int<N, T> r{};

#pragma unroll
  for (auto i = 0U; i < N; ++i) {
    auto aa = a[i];
    auto sum = aa + b[i];
    auto res = sum + carry;
    carry = (sum < aa) | (res < sum);
    r[i] = res;
  }

  auto reduced = subtract(r, modulus);
  T r_geq_modulus = reduced[N] ? 0 : 1;
  big_int<N, T> res = (carry + r_geq_modulus != 0) ? detail::first<N>(reduced) : r;
  return res;
}

template <typename T, size_t N>
constexpr auto mod_sub(big_int<N, T> a, big_int<N, T> b,
                       big_int<N, T> modulus) {
#ifdef __clang__
  if constexpr (std::is_same_v<T, uint64_t> && N == 4) {
    using u256 = _ExtInt(256);
    using u256_ar = big_int<N,T>;
    uint64_t carryin=0, carryout;
    uint64_t const* x = &a[0];
    uint64_t const* y = &b[0];
    u256_ar z;
    z[0] = __builtin_subcl(x[0], y[0], carryin, &carryout);
    carryin = carryout;
    z[1] = __builtin_subcl(x[1], y[1], carryin, &carryout);
    carryin = carryout;
    z[2] = __builtin_subcl(x[2], y[2], carryin, &carryout);
    carryin = carryout;
    z[3] = __builtin_subcl(x[3], y[3], carryin, &carryout);

    const uint64_t mask = -carryout;
    const u256_ar m256{mask, mask, mask, mask};
    const u256 pm = std::bit_cast<u256>(modulus) & std::bit_cast<u256>(m256);
    const auto retmp = std::bit_cast<u256_ar>(std::bit_cast<u256>(z)+pm);
    return retmp;
  }
#endif

  T carry{};
  big_int<N, T> r{};

  for (auto i = 0U; i < N; ++i) {
    auto aa = a[i];
    auto diff = aa - b[i];
    auto res = diff - carry;
    carry = (diff > aa) | (res > diff);
    r[i] = res;
  }

  auto adjusted_r = add_ignore_carry(r, modulus);
  big_int<N, T> res = carry ? adjusted_r : r;
  return res;
}



template <typename T, size_t N, T... Modulus>
constexpr auto mod_add(big_int<N, T> a, big_int<N, T> b, std::integer_sequence<T, Modulus...>) {
  big_int<sizeof...(Modulus), T> modulus{{Modulus...}};
  return mod_add(a, b, modulus);
}

template <typename T, size_t N1, size_t N2>
constexpr auto operator+(big_int<N1, T> a, big_int<N2, T> b) {
  return add(a, b);
}

template <typename T, size_t N1, size_t N2>
constexpr auto operator-(big_int<N1, T> a, big_int<N2, T> b) {
  return subtract(a, b);
}

template <typename T, size_t N, T... Modulus>
constexpr auto mod_sub(big_int<N, T> a, big_int<N, T> b, std::integer_sequence<T, Modulus...>) {
  big_int<sizeof...(Modulus), T> modulus{{Modulus...}};
  return mod_sub(a, b, modulus);
}



} // end namespace cbn

#endif
