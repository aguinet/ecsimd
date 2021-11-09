#ifndef ECSIMD_MODULAR_H
#define ECSIMD_MODULAR_H

#include <ecsimd/add.h>
#include <ecsimd/sub.h>

namespace ecsimd {

template <concepts::wide_bignum WBN>
auto mod_add(WBN const& a, WBN const& b, WBN const& p)
{
  const auto [sum, carry_add] = add(a,b);
  return sub_if_above(sum, p, !carry_add);
}

template <concepts::wide_bignum WBN>
auto mod_sub(WBN const& a, WBN const& b, WBN const& p)
{
  const auto [sum, carry_sub] = sub(a,b);
  return sub_if_above(sum, p, carry_sub);
}

} // ecsimd

#endif
