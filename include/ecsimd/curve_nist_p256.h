#ifndef ECSIMD_CURVE_NIST_P256_H
#define ECSIMD_CURVE_NIST_P256_H

#include <ecsimd/bignum.h>
#include <ecsimd/curve.h>
#include <ecsimd/serialization.h>
#include <ecsimd/literals.h>

using namespace ecsimd::literals;

namespace ecsimd {

// From https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-186-draft.pdf
struct curve_nist_p256 {
  using bn_type = bignum_256;
  struct P {
    static constexpr auto value = bn_from_bytes_BE<bn_type>("ffffffff00000001000000000000000000000000ffffffffffffffffffffffff"_hex);
  };
  struct A {
    // -3 % p
    static constexpr auto value = bn_from_bytes_BE<bn_type>("ffffffff00000001000000000000000000000000fffffffffffffffffffffffc"_hex);
  };
  struct B {
    static constexpr auto value = bn_from_bytes_BE<bn_type>("5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b"_hex);
  };
  struct Gx {
    static constexpr auto value = bn_from_bytes_BE<bn_type>("6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296"_hex);
  };
  struct Gy {
    static constexpr auto value = bn_from_bytes_BE<bn_type>("4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5"_hex);
  };
};

static_assert(concepts::wst_curve_am3<curve_nist_p256>);

} // ecsimd

#endif
