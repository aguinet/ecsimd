//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Project Contributors
  SPDX-License-Identifier: BSL-1.0
*/
//==================================================================================================
#pragma once

#include <eve/concept/compatible.hpp>
#include <eve/concept/value.hpp>
#include <eve/detail/apply_over.hpp>
#include <eve/detail/implementation.hpp>
#include <eve/detail/skeleton_calls.hpp>
#include <eve/module/core/decorator/numeric.hpp>
#include <eve/module/core/numeric/fma.hpp>
#include <eve/module/core/regular/all.hpp>

namespace eve::detail
{
template<real_value T, real_value U, real_value V>
EVE_FORCEINLINE auto
fms_(EVE_SUPPORTS(cpu_), numeric_type const&, T const& a, U const& b, V const& c) noexcept requires
    compatible_values<T, U> && compatible_values<T, V>
{
  return arithmetic_call(numeric(fms), a, b, c);
}

template<real_value T>
EVE_FORCEINLINE T
fms_(EVE_SUPPORTS(cpu_), numeric_type const&, T const& a, T const& b, T const& c) noexcept requires
    has_native_abi_v<T>
{
  return numeric(fma)(a, b, T(-c));
}
}