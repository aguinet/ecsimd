#ifndef FP_SIMD_BIG_INT
#define FP_SIMD_BIG_INT

#include <eve/wide.hpp>
#include <eve/product_type.hpp>

#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace ecsimd {

namespace details {

template <class BN, class LimbType, class Indices>
struct eve_struct_nlimbs_impl; 

template <class BN, class LimbType, size_t... I>
struct eve_struct_nlimbs_impl<BN, LimbType, std::integer_sequence<size_t, I...>> {

  template <size_t I_>
  using limb_type_ = LimbType;

  using type = eve::struct_support<BN,
        limb_type_<I>...>;
}; 

template <class BN, class LimbType, size_t NLimbs>
using eve_struct_nlimbs = typename eve_struct_nlimbs_impl<BN, LimbType,
      std::make_integer_sequence<size_t, NLimbs>>::type;

} // details

template <class LimbType, size_t NLimbs>
struct bignum: details::eve_struct_nlimbs<
                 bignum<LimbType, NLimbs>,
                 LimbType, NLimbs>
{
  static_assert(std::is_unsigned_v<LimbType>);
  using limb_type = LimbType;
  static constexpr size_t nlimbs = NLimbs;
};

using bignum_128 = bignum<uint32_t, 4>;
using bignum_256 = bignum<uint32_t, 8>;

} // ecsimd

#endif
