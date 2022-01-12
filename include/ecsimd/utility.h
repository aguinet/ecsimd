#ifndef ECSIMD_UTILITY_H
#define ECSIMD_UTILITY_H

#include <ecsimd/bignum.h>

#include <eve/function/store.hpp>
#include <eve/detail/meta.hpp>

#include <utility>
#include <type_traits>
#include <array>

namespace ecsimd {

namespace concepts {
template <class T>
concept array = std::is_same_v<T, std::array<typename T::value_type, std::tuple_size<T>::value>>;
} // concepts

namespace details {
// clang 13 does not support this concept-based constraints
template </*concepts::array */auto const& ar, class T, size_t... I>
consteval auto ar_to_is_impl(std::integer_sequence<size_t, I...>) {
    return std::integer_sequence<T, ar[I]...>{};
}
template </*concepts::array */auto const& ar>
consteval auto ar_to_is() {
    return ar_to_is_impl<ar,
      typename std::decay_t<decltype(ar)>::value_type>(std::make_integer_sequence<size_t, ar.size()>());
}
} // details

template </*concepts::array */auto const& ar>
using array_to_integer_sequence_t = decltype(details::ar_to_is<ar>());

template <std::unsigned_integral T, class C>
static auto wide_uasr(eve::wide<T, C> uv, auto S) {
  using Signed = eve::detail::make_integer_t<sizeof(T), signed>;
  using SW = eve::wide<Signed, C>;
  const auto sv = std::bit_cast<SW>(uv) >> S;
  return std::bit_cast<eve::wide<T,C>>(sv);
}

template <std::unsigned_integral T, class C>
static auto wide_mask_bit(eve::wide<T, C> v, auto B) {
  using wide_ty = eve::wide<T,C>;
  constexpr auto nbits = std::numeric_limits<T>::digits;

  const auto wlsb = wide_ty{1};
  return eve::logical<wide_ty>{wide_uasr(v << (nbits-B-1), nbits-1)};
}

} // ecsimd

#endif
