//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

#include <eve/constant/zero.hpp>
#include <eve/function/derivative.hpp>
#include <eve/function/floor.hpp>

namespace eve::detail
{

  template<floating_real_value T>
  EVE_FORCEINLINE constexpr T floor_(EVE_SUPPORTS(cpu_)
                                   , diff_type<1> const &
                                   , T x) noexcept
  {
    return zero(as(x));
  }
}
