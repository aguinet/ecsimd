#ifndef FP_SIMD_BIG_INT
#define FP_SIMD_BIG_INT

#include <eve/wide.hpp>
#include <eve/product_type.hpp>

#include <ctbignum/bigint.hpp>

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <concepts>
#include <bit>

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

template <std::unsigned_integral LimbType, size_t NLimbs>
struct bignum: details::eve_struct_nlimbs<
                 bignum<LimbType, NLimbs>,
                 LimbType, NLimbs>
{
  using limb_type = LimbType;
  using base_type = details::eve_struct_nlimbs<
                 bignum<LimbType, NLimbs>,
                 LimbType, NLimbs>;
  using cbn_type = cbn::big_int<NLimbs, LimbType>;

  static constexpr size_t nlimbs = NLimbs;

  static bignum from(limb_type const v0) {
    bignum ret;
    get<0>(ret) = v0;
    eve::detail::for_<1,1,nlimbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
      constexpr auto i = decltype(i_)::value;
      get<i>(ret) = 0;
    });
    return ret;
  }

  constexpr auto cbn() const {
    return std::bit_cast<cbn_type>(*this);
  }

  static constexpr bignum from(cbn_type const& v) {
    return std::bit_cast<bignum>(v);
  }

  template <eve::like<bignum> T>
  friend auto& operator&=(T& a, T const& b) {
    eve::detail::for_<0,1,nlimbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
      constexpr auto i = decltype(i_)::value;
      get<i>(a) &= get<i>(b);
    });
    return a;
  }

  template <eve::like<bignum> T>
  friend auto& operator|=(T& a, T const& b) {
    eve::detail::for_<0,1,nlimbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
      constexpr auto i = decltype(i_)::value;
      get<i>(a) |= get<i>(b);
    });
    return a;
  }

  template <eve::like<bignum> T>
  friend auto& operator^=(T& a, T const& b) {
    eve::detail::for_<0,1,nlimbs>([&](auto i_) EVE_LAMBDA_FORCEINLINE {
      constexpr auto i = decltype(i_)::value;
      get<i>(a) ^= get<i>(b);
    });
    return a;
  }
};

using bignum_128 = bignum<uint64_t, 2>;
using bignum_256 = bignum<uint64_t, 4>;
using bignum_512 = bignum<uint64_t, 8>;

template <class Bignum>
using wide_bignum = eve::wide<Bignum, eve::fixed<4>>;

namespace concepts {
template <class T>
concept bignum = std::same_as<T, bignum<typename T::limb_type, T::nlimbs>>;
template <class T>
concept wide_bignum = eve::simd_value<T> && bignum<typename T::value_type>;
template <class T>
concept bignum_cst = bignum<std::remove_cv_t<decltype(T::value)>>;
} // concepts

template <class>
struct bn_traits {};
template <concepts::bignum T>
struct bn_traits<T> {
  using bignum_type = T;
  using limb_type = typename T::limb_type;
  static constexpr size_t nlimbs = T::nlimbs;
};
template <concepts::wide_bignum T>
struct bn_traits<T>: public bn_traits<typename T::value_type> { };
template <concepts::bignum_cst P>
struct bn_traits<P>: public bn_traits<std::remove_cv_t<decltype(P::value)>> { };

template <class T>
using bn_t = typename bn_traits<T>::bignum_type;

template <class T>
using bn_limb_t = typename bn_traits<T>::limb_type;

template <class T>
static constexpr size_t bn_nlimbs = bn_traits<T>::nlimbs;

// Used as a return type of comparaisons and carries
template <concepts::wide_bignum WBN>
using cmp_res_t = eve::logical<eve::wide<bn_limb_t<WBN>, eve::cardinal_t<WBN>>>;

namespace concepts {
template <class T, class WBN>
concept cmp_res = std::same_as<T, cmp_res_t<WBN>>;
} // concepts

namespace details {
template <std::unsigned_integral LimbType, concepts::bignum BN>
constexpr auto remap_limb(BN const& bn) {
  using org_limb_type = bn_limb_t<BN>;
  constexpr auto org_nlimbs = bn_nlimbs<BN>;
  constexpr auto new_nlimbs = (sizeof(org_limb_type)*org_nlimbs)/sizeof(LimbType);
  // TODO: what if big-endian?
  return std::bit_cast<bignum<LimbType, new_nlimbs>>(bn);
}
} // details

template <concepts::bignum_cst P, std::unsigned_integral LimbType>
struct remap_limb
{
  struct type {
    static constexpr auto value = details::remap_limb<LimbType>(P::value);
  };
};

template <concepts::bignum_cst P, std::unsigned_integral LimbType>
using remap_limb_t = typename remap_limb<P, LimbType>::type;

template <concepts::wide_bignum WBN>
using wbn_zext_t = wide_bignum<bignum<bn_limb_t<WBN>, bn_nlimbs<WBN>*2>>;

} // ecsimd

#endif
