//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

#include <eve/constant/one.hpp>
#include <eve/function/derivative.hpp>
#include <eve/function/fnma.hpp>
#include <eve/function/minus.hpp>

namespace eve::detail
{
  template<floating_value T, floating_value U, floating_value V, auto N>
  EVE_FORCEINLINE  auto fnma_(EVE_SUPPORTS(cpu_)
                            , diff_type<N> const &
                            , T const &a
                            , U const &b
                            , V const &c) noexcept
  requires compatible_values<T, U>&&compatible_values<T, V>
  {
    return arithmetic_call(diff_type<N>{}(fnma), a, b, c);
  }

  template<floating_value T, auto N>
  EVE_FORCEINLINE  auto fnma_(EVE_SUPPORTS(cpu_)
                            , diff_type<N> const &
                            , T const &a
                            , T const &b
                            , T const &c) noexcept
  requires(has_native_abi_v<T>)
  {
    if constexpr(N == 1) return minus(b);
    else if constexpr(N == 2) return minus(a);
    else if constexpr(N == 3) return one(as(c));
  }
}
