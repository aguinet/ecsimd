#ifndef ECSIMD_CURVE_GROUPE_H
#define ECSIMD_CURVE_GROUPE_H

#include <ecsimd/curve.h>
#include <ecsimd/mgry_ops.h>
#include <ecsimd/gfp.h>
#include <ecsimd/jacobian_curve_point.h>
#include <ecsimd/swap.h>
#include <ecsimd/ifelse.h>

#include <optional>
#include <cassert>

#include <iostream>

namespace ecsimd {

template <class>
struct curve_group;

template <concepts::wst_curve_am3 Curve>
struct curve_group<Curve> {
  using WBN  = curve_wide_bn_t<Curve>;
  using BN = typename WBN::value_type;
  using WMBN = curve_wide_mgry_bn_t<Curve>;
  using gfp = GFp<typename Curve::P>;

  using WCP  = wide_curve_point<Curve>;
  using WJCP = wide_jacobian_curve_point<Curve>;

  static constexpr auto Bm = to_mgry<typename Curve::P>(Curve::B::value);
  static constexpr auto Am = to_mgry<typename Curve::P>(Curve::A::value);
  //static const WMBN Bwm;

  static auto WG() {
    return WCP{WBN{Curve::Gx::value}, WBN{Curve::Gy::value}};
  }

  static auto WJG() {
    return WJCP::from_affine(WG());
  }

  static std::optional<WMBN> compute_y(WMBN const& x) {
    // y^2 = x^3 + ax + b
    // a == -3
    const auto xpow3 = x.sqr() * x;
    const auto x3 = mgry_shift_left<1>(x) + x;
    const auto ypow2 = xpow3 + WMBN{WBN{Bm}} - x3;
    return ypow2.sqrt();
  }

  static std::optional<WBN> compute_y(WBN const& x) {
    const auto ret = compute_y(WMBN::from_classical(x));
    if (!ret) {
      return {};
    }
    return {ret->to_classical()};
  }

  // All the following co-Z Jacobian curve point computations are based on
  // https://eprint.iacr.org/2010/309.pdf.

  // Double-with-update (co-Z). P.z must be equal to mgry(1).
  static WJCP DBLU(WJCP& P) {
    // TODO: assert z == mgry(1)
    auto& X1 = P.x();
    auto& Y1 = P.y();

    const auto B = X1.sqr();
    const auto E = Y1.sqr();
    const auto L = E.sqr();
    const auto S = mgry_shift_left<1>((X1 + E).sqr() - B - L);
    const auto M = (mgry_shift_left<1>(B)+B) + WMBN{WBN{Am}};

    WJCP ret;
    ret.x() = M.sqr() - mgry_shift_left<1>(S);
    const auto Lm8 = mgry_shift_left<3>(L);
    ret.y() = M*(S-ret.x()) - Lm8;
    ret.z() = mgry_shift_left<1>(Y1);

    // Update P
    X1 = S;
    Y1 = Lm8;
    P.z() = ret.z();

    return ret;
  }

  // Co-Z Addition with updates. Return P + O, and updates P to get the same
  // Z coordinate as the returned point.
  static WJCP ZADDU(WJCP& P, WJCP const& O) {
    assert(eve::all(P.z().wbn() == O.z().wbn()));
    auto& Z = P.z();
    auto& X1 = P.x();
    auto& Y1 = P.y();
    auto const& X2 = O.x();
    auto const& Y2 = O.y();

    const auto C = (X1-X2).sqr();
    const auto W1 = X1*C;
    const auto W2 = X2*C;
    const auto D = (Y1-Y2).sqr();
    const auto A1 = Y1*(W1-W2);

    WJCP ret;
    ret.x() = D-W1-W2;
    ret.y() = (Y1-Y2)*(W1-ret.x())-A1;
    ret.z() = Z*(X1-X2);

    // Update P
    X1 = W1;
    Y1 = A1;
    Z  = ret.z();

    return ret;
  }

  // Co-Z doubling-addition with updates. Returns (2P+Q) and updates Q with
  // same Z than returned point. 
  static WJCP ZDAU(WJCP const& P, WJCP& Q) {
    assert(eve::all(P.z().wbn() == Q.z().wbn()));
    const auto& X1 = P.x();
    const auto& Y1 = P.y();
    const auto& Z  = P.z();

    auto& X2 = Q.x();
    auto& Y2 = Q.y();

    const auto Cp = (X1 - X2).sqr();
    const auto W1p = X1*Cp;
    const auto W2p = X2*Cp;
    const auto Dp = (Y1-Y2).sqr();
    const auto A1p = Y1*(W1p-W2p);
    const auto X3pc = Dp - W1p - W2p;
    const auto C = (X3pc - W1p).sqr();
    const auto Y3p = ((Y1-Y2) + (W1p - X3pc)).sqr() - Dp - C - mgry_shift_left<1>(A1p);
    const auto W1 = mgry_shift_left<2>(X3pc)*C;
    const auto W2 = mgry_shift_left<2>(W1p)*C;
    const auto D = (Y3p-mgry_shift_left<1>(A1p)).sqr();
    const auto A1 = Y3p*(W1-W2);

    WJCP ret;
    ret.x() = D-W1-W2;
    ret.y() = (Y3p-mgry_shift_left<1>(A1p))*(W1-ret.x()) - A1;
    ret.z() = Z*((X1-X2+X3pc-W1p).sqr() - Cp - C);

    const auto Dc = (Y3p+mgry_shift_left<1>(A1p)).sqr();
    X2 = Dc - W1 - W2;
    Y2 = (Y3p + mgry_shift_left<1>(A1p))*(W1-X2)-A1;
    Q.z() = ret.z();

    return ret;
  }

  static WJCP ADD_Z2_1(WJCP const& A, WJCP const& B) {
    // TODO: assert B.z() == R
    const auto& X1 = A.x();
    const auto& Y1 = A.y();
    const auto& Z1 = A.z();
    const auto& X2 = B.x();
    const auto& Y2 = B.y();

    const auto Z1Z1 = Z1.sqr();
    const auto U2 = X2*Z1Z1;
    const auto S2 = Y2*Z1*Z1Z1;
    const auto H = U2-X1;
    const auto HH = H.sqr();
    const auto I = mgry_shift_left<2>(HH);
    const auto J = H*I;
    const auto r = mgry_shift_left<1>(S2-Y1);
    const auto V = X1*I;

    WJCP ret;
    ret.x() = r.sqr()-J-mgry_shift_left<1>(V);
    ret.y() = r*(V-ret.x())-mgry_shift_left<1>(Y1)*J;
    ret.z() = (Z1+H).sqr()-Z1Z1-HH;

    return ret;
  }

  // Co-Z point triple. Returns 3*P, and updates P to get the same Z coordinate
  // as the returned point.
  static WJCP TRPLU(WJCP& P) {
    const auto dbl = DBLU(P);
    return ZADDU(P, dbl);
  }

  // Scalar multiplication, 4 scalars x 4 points
  static WJCP scalar_mult(WBN const& x, WJCP P) {
    using limb_t = bn_limb_t<WBN>;
    using wide_limb_t = eve::wide<limb_t, eve::cardinal_t<WBN>>;
    constexpr auto limb_nbits = std::numeric_limits<limb_t>::digits;

    const auto oppP = P.opposite();
    auto base = TRPLU(P);
    swap_if_same_z(wide_mask_bit(eve::get<0>(x), 1), P, base);

    eve::detail::for_<0,1,bn_nlimbs<WBN>>([&](auto l_) ECSIMD_LAMBDA_FORCEINLINE {
      constexpr auto l = decltype(l_)::value;
      auto const& lx = eve::get<l>(x);

      size_t b = 0;
      if constexpr (l == 0) {
        b = 2;
      }
      for (; b < limb_nbits; ++b) {
        const auto mask = wide_mask_bit(lx, b);
        swap_if_same_z(mask, P, base);
        base = ZDAU(base, P);
        swap_if_same_z(mask, P, base);
      }
    });

    const auto wlsb = wide_limb_t{1};
    const auto meven = (eve::get<0>(x) & wlsb) == eve::zero(eve::as<wide_limb_t>());
    auto Psub = ADD_Z2_1(P, oppP);
    return if_else(meven, Psub, P);
  }

  // Scalar multiplication, 1 scalar x 4 points
  static WJCP scalar_mult_1s(BN const& x, WJCP P) {
    using limb_t = bn_limb_t<WBN>;
    using wide_limb_t = eve::wide<limb_t, eve::cardinal_t<WBN>>;
    constexpr auto limb_nbits = std::numeric_limits<limb_t>::digits;

    const auto oppP = P.opposite();
    auto base = TRPLU(P);

    auto xbit = (kumi::get<0>(x) >> 1) & 1;
    WJCP* R[2];
    R[1-xbit] = &base;
    R[xbit] = &P;

    eve::detail::for_<0,1,bn_nlimbs<WBN>>([&](auto l_) ECSIMD_LAMBDA_FORCEINLINE {
      constexpr auto l = decltype(l_)::value;
      auto const lx = kumi::get<l>(x);

      size_t b = 0;
      if constexpr (l == 0) {
        b = 2;
      }
      for (; b < limb_nbits; ++b) {
        xbit = (lx >> b) & 1;
        *R[1-xbit] = ZDAU(*R[1-xbit], *R[xbit]);
      }
    });

    xbit = kumi::get<0>(x) & 1;
    *R[xbit] = ADD_Z2_1(*R[0], oppP);
    return *R[0];
  }
};

//template <concepts::wst_curve_am3 Curve>
//const curve_wide_mgry_bn_t<Curve> curve_group<Curve>::Bwm = curve_wide_mgry_bn_t<Curve>{curve_wide_bn_t<Curve>{curve_group<Curve>::Bm}};

} // ecsimd

#endif
