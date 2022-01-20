#include <cstdint>
#include <bit>
#include <limits>
#include <cassert>

#include <ctbignum/bigint.hpp>
#include <ctbignum/decimal_literals.hpp>
#include <ctbignum/montgomery.hpp>
#include <ctbignum/io.hpp>
#include <ctbignum/mod_inv.hpp>
#include <ctbignum/invariant_div.hpp>

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

#if 0
__attribute__((noinline)) cbn_u256 fast_reduce(cbn_u512 const& n) {
  //return cbn::mod(n, P);
  using cbn_u512_l32 = cbn::big_int<16, uint32_t>;
  const auto n32 = std::bit_cast<cbn_u512_l32>(n);
  cbn::big_int<16, uint64_t> a;
  std::transform(n32.begin(), n32.end(), a.begin(), [](uint32_t const v) -> uint64_t { return static_cast<uint64_t>(v); });

  const uint64_t l0 =              a[0]+a[8]+a[9]+2*a[10]+3*a[11]+5*a[12]+9*a[13]+15*a[14]+26*a[15];
  const uint64_t l1 = (l0 >> 32) + a[1]+a[9]+a[10]+2*a[11]+3*a[12]+5*a[13]+9*a[14]+15*a[15];
  const uint64_t l2 = (l1 >> 32) + a[2]+a[10]+a[11]+2*a[12]+3*a[13]+5*a[14]+9*a[15];
  const uint64_t l3 = (l2 >> 32) + a[3]+a[8]+a[9]+2*a[10]+4*a[11]+6*a[12]+11*a[13]+18*a[14]+31*a[15];
  const uint64_t l4 = (l3 >> 32) + a[4]+a[9]+a[10]+2*a[11]+4*a[12]+6*a[13]+11*a[14]+18*a[15];
  const uint64_t l5 = (l4 >> 32) + a[5]+a[10]+a[11]+2*a[12]+4*a[13]+6*a[14]+11*a[15];
  const uint64_t l6 = (l5 >> 32) + a[6]+a[8]+a[9]+2*a[10]+4*a[11]+6*a[12]+11*a[13]+19*a[14]+32*a[15];
  const uint64_t l7 = (l6 >> 32) + a[7]+a[8]+2*a[9]+3*a[10]+5*a[11]+9*a[12]+15*a[13]+26*a[14]+45*a[15];

  const uint64_t r0 = (l0&0xffffffff) | (l1<<32);
  const uint64_t r1 = (l2&0xffffffff) | (l3<<32);
  const uint64_t r2 = (l4&0xffffffff) | (l5<<32);
  const uint64_t r3 = (l5&0xffffffff) | (l7<<32);
  const uint64_t r4 = l7>>32;

  // Subtract P
  cbn::big_int<5,uint64_t> ret{r0,r1,r2,r3,r4};
  return cbn::detail::first<4>(ret);
  //return cbn::detail::first<4>(cbn::mod(ret, P));
}
#endif

cbn_u256 fast_reduce(cbn_u512 const& n);
cbn_u256 mul_mod(cbn_u256 const& a, cbn_u256 const& b) {
  const auto m = cbn::mul(a,b);
  return fast_reduce(m);
}

cbn_u256 sqr_mod(cbn_u256 const& a) {
  return mul_mod(a,a);
}

struct gfp {
  gfp() = default;

  constexpr gfp(cbn_u256 const& n):
    n_(n)
  { }

  gfp sqr() const {
    return gfp{sqr_mod(n_)};
  }

  gfp opposite() const {
    assert(n_ < cbn_P);
    return gfp{cbn::subtract_ignore_carry(cbn_P,n_)};
  }

  cbn_u256 const& bn() const { return n_; }
  cbn_u256& bn() { return n_; }

  auto operator<=>(gfp const& o) const {
    return bn() <=> o.bn();
  }

private:
  cbn_u256 n_;
};

gfp operator*(gfp const& a, gfp const& b) {
  return gfp{mul_mod(a.bn(), b.bn())};
}

gfp operator+(gfp const& a, gfp const& b) {
  return gfp{cbn::mod_add(a.bn(), b.bn(), cbn_P)};
}

gfp operator-(gfp const& a, gfp const& b) {
  return gfp{cbn::mod_sub(a.bn(), b.bn(), cbn_P)};
}

template <size_t Count>
gfp gfp_shift_left(gfp const& a) {
  // TODO: optimize
  gfp ret = a;
#pragma unroll
  for (size_t i = 0; i < Count; ++i) {
    ret = ret+ret;
  }
  return ret;
}

cbn_u256 fast_reduce(cbn_u512 const& n) {
  using cbn_u512_l32 = cbn::big_int<16, uint32_t>;
  using cbn_u256_l32 = cbn::big_int<8, uint32_t>;
  const auto a = std::bit_cast<cbn_u512_l32>(n);

  const gfp s2(std::bit_cast<cbn_u256>(cbn_u256_l32{0,0,0,a[11],a[12],a[13],a[14],a[15]}));
  const gfp s3(std::bit_cast<cbn_u256>(cbn_u256_l32{0,0,0,a[12],a[13],a[14],a[15],0}));
  const gfp s4(std::bit_cast<cbn_u256>(cbn_u256_l32{a[8],a[9],a[10],0,0,0,a[14],a[15]}));
  const gfp s5(std::bit_cast<cbn_u256>(cbn_u256_l32{a[9],a[10],a[11],a[13],a[14],a[15],a[13],a[8]}));
  const gfp s6(std::bit_cast<cbn_u256>(cbn_u256_l32{a[11],a[12],a[13],0,0,0,a[8],a[10]}));
  const gfp s7(std::bit_cast<cbn_u256>(cbn_u256_l32{a[12],a[13],a[14],a[15],0,0,a[9],a[11]}));
  const gfp s8(std::bit_cast<cbn_u256>(cbn_u256_l32{a[14],a[14],a[15],a[8],a[9],a[10],0,a[12]}));
  const gfp s9(std::bit_cast<cbn_u256>(cbn_u256_l32{a[14],a[15],0,a[9],a[10],a[11],0,a[13]}));

  return (cbn::detail::first<4>(n)+gfp_shift_left<1>(s2)+gfp_shift_left<1>(s3)+s4+s5-s6-s7-s8-s9).bn();
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
    ret.x_ = gfp{pt.x()};
    ret.y_ = gfp{pt.y()};
    ret.z_ = gfp{cbn_u256{1,0,0,0}};
    return ret;
  }

  curve_point to_affine() const {
    // Unoptimized, only for testing purposes
    const auto& Z = z_.bn();
    const auto invZ = cbn::mod_inv(Z, cbn_P);
    const auto invZ2 = sqr_mod(invZ);
    const auto invZ3 = mul_mod(invZ2,invZ);

    curve_point ret;
    ret.x() = mul_mod(x_.bn(), invZ2);
    ret.y() = mul_mod(y_.bn(), invZ3);
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
  gfp x_;
  gfp y_;
  gfp z_;
};

using JCP = jacobian_curve_point;

// -3 % p
constexpr gfp Am = gfp{cbn::to_big_int(115792089210356248762697446949407573530086143415290314195533631308867097853948_Z)};

// Double-with-update (co-Z). P.z must be equal to mgry(1).
JCP DBLU(JCP& P) {
  // TODO: assert z == mgry(1)
  auto& X1 = P.x();
  auto& Y1 = P.y();

  const auto B = X1.sqr();
  const auto E = Y1.sqr();
  const auto L = E.sqr();
  const auto S = gfp_shift_left<1>((X1 + E).sqr() - B - L);
  const auto M = (gfp_shift_left<1>(B)+B) + Am;

  JCP ret;
  ret.x() = M.sqr() - gfp_shift_left<1>(S);
  const auto Lm8 = gfp_shift_left<3>(L);
  ret.y() = M*(S-ret.x()) - Lm8;
  ret.z() = gfp_shift_left<1>(Y1);

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
  const auto Y3p = ((Y1-Y2) + (W1p - X3pc)).sqr() - Dp - C - gfp_shift_left<1>(A1p);
  const auto W1 = gfp_shift_left<2>(X3pc)*C;
  const auto W2 = gfp_shift_left<2>(W1p)*C;
  const auto D = (Y3p-gfp_shift_left<1>(A1p)).sqr();
  const auto A1 = Y3p*(W1-W2);

  JCP ret;
  ret.x() = D-W1-W2;
  ret.y() = (Y3p-gfp_shift_left<1>(A1p))*(W1-ret.x()) - A1;
  ret.z() = Z*((X1-X2+X3pc-W1p).sqr() - Cp - C);

  const auto Dc = (Y3p+gfp_shift_left<1>(A1p)).sqr();
  X2 = Dc - W1 - W2;
  Y2 = (Y3p + gfp_shift_left<1>(A1p))*(W1-X2)-A1;
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
  const auto I = gfp_shift_left<2>(HH);
  const auto J = H*I;
  const auto r = gfp_shift_left<1>(S2-Y1);
  const auto V = X1*I;

  JCP ret;
  ret.x() = r.sqr()-J-gfp_shift_left<1>(V);
  ret.y() = r*(V-ret.x())-gfp_shift_left<1>(Y1)*J;
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

template <class T>
T rnd_bn() {
  T ret;
  getentropy(&ret, sizeof(ret));
  return ret;
}

cbn_u256 rnd_u256() {
  return rnd_bn<cbn_u256>();
}

cbn_u512 rnd_u512() {
  return rnd_bn<cbn_u512>();
}

void bench_redc(benchmark::State& S) {
  const cbn_u256 a = cbn::mod(rnd_u256(), P);
  const cbn_u256 b = cbn::mod(rnd_u256(), P);
  const cbn_u512 m = cbn::mul(a,b);

  for (auto _: S) {
    const auto red = fast_reduce(m);
    benchmark::DoNotOptimize(red);
  }
}

void bench_p256(benchmark::State& S) {
  const cbn_u256 x = cbn::mod(rnd_u256(), P);

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
  const cbn_u256 Gx = cbn::to_big_int(48439561293906451759052585252797914202762949526041747995844080717082404635286_Z);
  const cbn_u256 Gy = cbn::to_big_int(36134250956749795798585127919587881956611106672985015071877198253568414405109_Z);
  const curve_point G(Gx,Gy);
  const auto JG = jacobian_curve_point::from_affine(G);

  const auto a = cbn::mod(rnd_u256(), P);
  const auto b = cbn::mod(rnd_u256(), P);
  const auto m = a*b;
  std::cout << "m= " << m << std::endl;
  const auto redref = cbn::mod(m,P);
  const auto redfast = fast_reduce(m);
  std::cout << "m%p        = " << std::hex << redref << std::dec << std::endl;
  std::cout << "fast_red(m)= " << std::hex << redfast << std::dec << std::endl;
  std::cout << (redref == redfast) << std::endl;

  const auto P = scalar_mult(a, JG).to_affine();
  std::cout << "a=" << a << std::endl;
  std::cout << "x=" << P.x() << std::endl;
  std::cout << "y=" << P.y() << std::endl;

  benchmark::RegisterBenchmark("redc", bench_redc);
  benchmark::RegisterBenchmark("scalar_mult_p256", bench_p256)->Unit(benchmark::kMicrosecond);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
