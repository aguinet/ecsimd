//==================================================================================================
/*
  EVE - Expressive Vector Engine
  Copyright : EVE Contributors & Maintainers
  SPDX-License-Identifier: MIT
*/
//==================================================================================================
#pragma once

#include <eve/function/diff/log.hpp>
#include <eve/function/derivative.hpp>
#include <eve/function/is_gtz.hpp>
#include <eve/function/rec.hpp>

namespace eve::detail
{

  template<floating_real_value T>
  EVE_FORCEINLINE constexpr T log10_(EVE_SUPPORTS(cpu_)
                                  , diff_type<1> const &
                                  , T const &x) noexcept
  {
    auto invlog10 = T(0.4342944819032518276511289189);
    return if_else(is_gtz(x), rec(x)*invlog10, allbits); ;
  }
}
