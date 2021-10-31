#ifndef FPSIMD_BIGINT_ADD_H
#define FPSIMD_BIGINT_ADD_H

#include <eve/wide.hpp>
#include <ecsimd/bignum.h>

namespace ecsimd {

template <class T, size_t N>
static auto add(eve::wide<bignum<T, N>> const& a, eve::wide<bignum<T, N>> const& b)
{
  using limb_type = T;
  constexpr auto nlimbs = N;
  eve::wide<limb_type> carry; 
  eve::wide<bignum<T,N>> ret;
  eve::detail::for_<0, 1, nlimbs>([&](auto i_) {
    constexpr auto i = decltype(i_)::value;
    const auto aa = eve::get<i>(a);
    const auto bb = eve::get<i>(b);
    const auto sum = aa + bb;
    if constexpr (i == 0) {
      carry = -(sum < aa).mask();
      eve::get<i>(ret) = sum;
    }
    else {
      const auto res = sum + carry;
      carry = -((sum < aa) || (res < sum)).mask();
      eve::get<i>(ret) = res;
    }
  });

  return ret;
}

} // ecsimd

#endif
