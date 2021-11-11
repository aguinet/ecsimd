#ifndef ECSIMD_UTILITY_H
#define ECSIMD_UTILITY_H

#include <ecsimd/bignum.h>

#include <eve/function/store.hpp>

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

} // ecsimd

#endif
