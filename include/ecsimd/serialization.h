#ifndef FPSIMD_SERIALIZATION_H
#define FPSIMD_SERIALIZATION_H

#include <ecsimd/bignum.h>
#include <ecsimd/intmem.h>

#include <cstdint>
#include <array>

namespace ecsimd {

template <class Bignum>
static constexpr Bignum bn_from_bytes_BE(uint8_t const* bytes)
{
  using limb_type = typename Bignum::limb_type;
  constexpr size_t nlimbs = Bignum::nlimbs;
  Bignum ret;
  kumi::for_each_index( [&](auto i_, auto& m) {
    constexpr auto i = decltype(i_)::value;
    m = intmem::loadu_be<limb_type>(&bytes[(nlimbs-i-1)*sizeof(limb_type)]);
  }, ret);
  return ret;
}

template <class Bignum>
static constexpr Bignum bn_from_bytes_BE(std::array<uint8_t, sizeof(Bignum)> const& bytes)
{
  return bn_from_bytes_BE<Bignum>(&bytes[0]);
}

template <class Bignum>
static void bn_to_bytes_BE(uint8_t* out, Bignum const& v)
{
  using limb_type = typename Bignum::limb_type;
  constexpr size_t nlimbs = Bignum::nlimbs;
  kumi::for_each_index( [&](auto i_, auto const& m) {
    constexpr auto i = decltype(i_)::value;
    intmem::storeu_be<limb_type>(&out[(nlimbs-i-1)*sizeof(limb_type)], m);
  }, v);
}

template <class Bignum>
static auto bn_to_bytes_BE(Bignum const& v)
{
  std::array<uint8_t, sizeof(Bignum)> Ret;
  bn_to_bytes_BE(&Ret[0], v);
  return Ret;
}

} // ecsimd

#endif
