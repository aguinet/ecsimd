#include <cstdint>
#include <bit>
#include <limits>
#include <cassert>

#include <ctbignum/bigint.hpp>
#include <ctbignum/decimal_literals.hpp>
#include <ctbignum/montgomery.hpp>
#include <ctbignum/io.hpp>
#include <ctbignum/mod_inv.hpp>

#include <benchmark/benchmark.h>

#include <unistd.h>
#include <iostream>

using namespace cbn::literals;

constexpr size_t nlimbs = 4;
using limb_type = uint64_t;

using u256 = _ExtInt(256);
using cbn_u256 = cbn::big_int<4, uint64_t>;

using u512 = _ExtInt(512);
using cbn_u512 = cbn::big_int<8, uint64_t>;

namespace {


constexpr u256 from_cbn(cbn_u256 const& v) {
  return std::bit_cast<u256>(v);
}

constexpr cbn_u256 to_cbn(u256 const& v) {
  return std::bit_cast<cbn_u256>(v);
}

constexpr cbn_u512 to_cbn_512(u512 const& v) {
  return std::bit_cast<cbn_u512>(v);
}

constexpr auto P = 115792089210356248762697446949407573530086143415290314195533631308867097853951_Z;
constexpr auto cbn_P = cbn::to_big_int(P);
constexpr auto R   = cbn::detail::unary_encoding<nlimbs, nlimbs + 1, limb_type>();
constexpr auto Rsq = cbn::detail::unary_encoding<2 * nlimbs, 2 * nlimbs + 1, limb_type>();
constexpr auto R_p   = cbn::div(R, P).remainder;
constexpr auto Rsq_p = cbn::div(Rsq, P).remainder;
constexpr auto Pm1_by_R_p = cbn::div((cbn_P-cbn_u256{1})*R, P).remainder;

constexpr cbn_u256 mgry_reduce(cbn_u512 const& n) {
  return cbn::montgomery_reduction(n, P);
}

constexpr cbn_u256 mgry_mul(cbn_u256 const& a, cbn_u256 const& b) {
  return mgry_reduce(cbn::mul(a,b));
}

cbn_u256 mgry_add(cbn_u256 const& a, cbn_u256 const& b) {
  return cbn::mod_add(a,b,cbn_P);
}

cbn_u256 mgry_sub(cbn_u256 const& a, cbn_u256 const& b) {
  return cbn::mod_sub(a,b,cbn_P);
}

struct mgry {
  mgry() = default;

  constexpr mgry(cbn_u256 const& n):
    n_(n)
  { }

  static constexpr mgry from_classical(cbn_u256 const& n) {
    return mgry{mgry_mul(n, Rsq_p)};
  }

  mgry sqr() const {
    return mgry_mul(n_,n_);
  }

  cbn_u256 to_classical() const {
    const cbn_u512 wn = cbn::detail::pad<nlimbs>(n_);
    return mgry_reduce(wn);
  }

  mgry opposite() const {
    auto ret = mgry_sub(n_, R_p);
    return mgry{mgry_sub(Pm1_by_R_p, ret)};
  }

  cbn_u256 const& bn() const { return n_; }
  cbn_u256& bn() { return n_; }

  auto operator<=>(mgry const& o) const {
    return bn() <=> o.bn();
  }

private:
  cbn_u256 n_;
};

mgry operator*(mgry const& a, mgry const& b) {
  return mgry{mgry_mul(a.bn(), b.bn())};
}

mgry operator+(mgry const& a, mgry const& b) {
  return mgry{mgry_add(a.bn(), b.bn())};
}

mgry operator-(mgry const& a, mgry const& b) {
  return mgry{mgry_sub(a.bn(), b.bn())};
}

template <size_t Count>
mgry mgry_shift_left(mgry const& a) {
  // TODO: optimize
  mgry ret = a;
#pragma unroll
  for (size_t i = 0; i < Count; ++i) {
    ret = ret+ret;
  }
  return ret;
}

struct curve_point {
  curve_point() = default;
  curve_point(curve_point const&) = default;

  curve_point(cbn_u256 const& x, cbn_u256 const& y):
    x_(x), y_(y)
  { }

  auto& x() { return x_; }
  auto& y() { return y_; }

  auto const& x() const { return x_; }
  auto const& y() const { return y_; }


private:
  cbn_u256 x_;
  cbn_u256 y_;
};

struct jacobian_curve_point {
  static jacobian_curve_point from_affine(curve_point const& pt) {
    jacobian_curve_point ret;
    ret.x_ = mgry::from_classical(pt.x());
    ret.y_ = mgry::from_classical(pt.y());
    ret.z_ = mgry{R_p};
    return ret;
  }

  curve_point to_affine() const {
    // Unoptimized, only for testing purposes
    const auto Z = z_.to_classical();
    const auto invZ = cbn::mod_inv(Z, cbn_P);
    const auto invZ2 = (invZ*invZ)%cbn_P;
    const auto invZ3 = (invZ2*invZ)%cbn_P;

    curve_point ret;
    ret.x() = (x_.to_classical() * invZ2)%cbn_P;
    ret.y() = (y_.to_classical() * invZ3)%cbn_P;
    return ret;
  }

  auto& x() { return x_; }
  auto& y() { return y_; }
  auto& z() { return z_; }

  auto const& x() const { return x_; }
  auto const& y() const { return y_; }
  auto const& z() const { return z_; }

  jacobian_curve_point opposite() const {
    jacobian_curve_point ret;
    ret.x() = x();
    ret.y() = y().opposite();
    ret.z() = z();
    return ret;
  }


private:
  mgry x_;
  mgry y_;
  mgry z_;
};

using JCP = jacobian_curve_point;

// mgry(-3, p)
constexpr mgry Am = mgry::from_classical(cbn::to_big_int(115792089210356248762697446949407573530086143415290314195533631308867097853948_Z));

// Double-with-update (co-Z). P.z must be equal to mgry(1).
JCP DBLU(JCP& P) {
  // TODO: assert z == mgry(1)
  auto& X1 = P.x();
  auto& Y1 = P.y();

  const auto B = X1.sqr();
  const auto E = Y1.sqr();
  const auto L = E.sqr();
  const auto S = mgry_shift_left<1>((X1 + E).sqr() - B - L);
  const auto M = (mgry_shift_left<1>(B)+B) + Am;

  JCP ret;
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
JCP ZADDU(JCP& P, JCP const& O) {
  assert(P.z().bn() == O.z().bn());
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

  JCP ret;
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
JCP ZDAU(JCP const& P, JCP& Q) {
  assert(P.z().bn() == O.z().bn());
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

  JCP ret;
  ret.x() = D-W1-W2;
  ret.y() = (Y3p-mgry_shift_left<1>(A1p))*(W1-ret.x()) - A1;
  ret.z() = Z*((X1-X2+X3pc-W1p).sqr() - Cp - C);

  const auto Dc = (Y3p+mgry_shift_left<1>(A1p)).sqr();
  X2 = Dc - W1 - W2;
  Y2 = (Y3p + mgry_shift_left<1>(A1p))*(W1-X2)-A1;
  Q.z() = ret.z();

  return ret;
}

JCP ADD_Z2_1(JCP const& A, JCP const& B) {
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

  JCP ret;
  ret.x() = r.sqr()-J-mgry_shift_left<1>(V);
  ret.y() = r*(V-ret.x())-mgry_shift_left<1>(Y1)*J;
  ret.z() = (Z1+H).sqr()-Z1Z1-HH;

  return ret;
}

// Co-Z point triple. Returns 3*P, and updates P to get the same Z coordinate
// as the returned point.
JCP TRPLU(JCP& P) {
  const auto dbl = DBLU(P);
  return ZADDU(P, dbl);
}

static JCP scalar_mult(cbn_u256 const& x, JCP P) {
  using limb_t = limb_type;
  constexpr auto limb_nbits = std::numeric_limits<limb_t>::digits;

  const auto oppP = P.opposite();
  auto base = TRPLU(P);

  auto xbit = (x[0] >> 1) & 1;
  JCP* R[2];
  R[1-xbit] = &base;
  R[xbit] = &P;

#pragma unroll
  for (size_t l = 0; l < nlimbs; ++l) {
    auto const lx = x[l];

    size_t b = 0;
    if (l == 0) {
      b = 2;
    }
    for (; b < limb_nbits; ++b) {
      xbit = (lx >> b) & 1;
      *R[1-xbit] = ZDAU(*R[1-xbit], *R[xbit]);
    }
  }

  xbit = x[0] & 1;
  *R[xbit] = ADD_Z2_1(*R[0], oppP);
  return *R[0];
}


} // anonymous

namespace {

cbn_u256 rnd_u256() {
cbn_u256 ret;
getentropy(&ret, sizeof(ret));
return ret;
}

void bench_p256(benchmark::State& S) {
  const cbn_u256 x = rnd_u256();

  const cbn_u256 Gx = cbn::to_big_int(48439561293906451759052585252797914202762949526041747995844080717082404635286_Z);
  const cbn_u256 Gy = cbn::to_big_int(36134250956749795798585127919587881956611106672985015071877198253568414405109_Z);
  const curve_point G(Gx,Gy);
  const auto JG = jacobian_curve_point::from_affine(G);

  for (auto _: S) {
    const auto P = scalar_mult(x, JG);
    benchmark::DoNotOptimize(P);
  }
}

} // namespace


int main(int argc, char** argv)
{
  const cbn_u256 a = rnd_u256();
  const cbn_u256 b = rnd_u256();

  const auto ma = mgry::from_classical(a);
  const auto mb = mgry::from_classical(b);

  const auto mref = cbn::mul(a,b)%cbn_P;
  const auto mmul = ma*mb;
  const cbn_u256 mul = mmul.to_classical();

  std::cout << (mref == mul) << std::endl;
  std::cout << ((ma+mb).to_classical() == cbn::mod_add(a,b,cbn_P)) << std::endl;
  std::cout << ((ma-mb).to_classical() == cbn::mod_sub(a,b,cbn_P)) << std::endl;

  const cbn_u256 Gx = cbn::to_big_int(48439561293906451759052585252797914202762949526041747995844080717082404635286_Z);
  const cbn_u256 Gy = cbn::to_big_int(36134250956749795798585127919587881956611106672985015071877198253568414405109_Z);
  const curve_point G(Gx,Gy);
  const auto JG = jacobian_curve_point::from_affine(G);

  const auto P = scalar_mult(a, JG).to_affine();
  std::cout << "a=" << a << std::endl;
  std::cout << "x=" << P.x() << std::endl;
  std::cout << "y=" << P.y() << std::endl;

  benchmark::RegisterBenchmark("scalar_mult_p256", bench_p256)->Unit(benchmark::kMicrosecond);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
