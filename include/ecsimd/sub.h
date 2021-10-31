#ifndef FPSIMD_BIGINT_SUB_H
#define FPSIMD_BIGINT_SUB_H

#include <eve/wide.hpp>
#include <ecsimd/bignum.h>

namespace ecsimd {

template <class T, size_t N>
static auto sub(eve::wide<bignum<T, N>> const& a, eve::wide<bignum<T, N>> const& b)
{
  using limb_type = T;
  constexpr auto nlimbs = N;
  eve::wide<limb_type> carry; 
  eve::wide<bignum<T,N>> ret;
  eve::detail::for_<0, 1, nlimbs>([&](auto i_) {
    constexpr auto i = decltype(i_)::value;
    const auto aa = eve::get<i>(a);
    const auto bb = eve::get<i>(b);
    const auto diff = aa - bb;
    if constexpr (i == 0) {
      carry = -(diff > aa).mask();
      eve::get<i>(ret) = diff;
    }
    else {
      const auto res = diff - carry;
      carry = -((diff > aa) || (res > diff)).mask();
      eve::get<i>(ret) = res;
    }
  });

  return ret;
}

} // ecsimd

#endif
