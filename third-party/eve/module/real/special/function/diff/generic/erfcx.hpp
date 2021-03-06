//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

#include <eve/function/erfcx.hpp>
#include <eve/function/derivative.hpp>
#include <eve/function/fma.hpp>

namespace eve::detail
{

  template<floating_real_value T>
  EVE_FORCEINLINE constexpr T erfcx_(EVE_SUPPORTS(cpu_)
                                  , diff_type<1> const &
                                  , T const &x) noexcept
  {
    auto twoosqrtpi = T(1.1283791670955125738961589);
    return fms(2*x, erfcx(x),  twoosqrtpi);
  }
}
