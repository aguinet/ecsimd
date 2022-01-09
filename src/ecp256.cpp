#include <cstdint>
#include <bit>
#include <limits>

#include <ctbignum/bigint.hpp>
#include <ctbignum/decimal_literals.hpp>
#include <ctbignum/montgomery.hpp>
#include <ctbignum/io.hpp>

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
constexpr auto R   = cbn::detail::unary_encoding<nlimbs, nlimbs + 1, limb_type>();
constexpr auto Rsq = cbn::detail::unary_encoding<2 * nlimbs, 2 * nlimbs + 1, limb_type>();
constexpr auto Rsq_p = cbn::div(Rsq, P).remainder;

static cbn_u256 mgry_reduce(cbn_u512 const& n) {
  return cbn::montgomery_reduction(n, P);
}

cbn_u256 mgry_mul(cbn_u256 const& a, cbn_u256 const& b) {
  return mgry_reduce(cbn::mul(a,b));
  const u256 wa = from_cbn(a);
  const u256 wb = from_cbn(b);
  return mgry_reduce(to_cbn_512(wa*wb));
}

struct mgry_u256 {
  mgry_u256() = default;

  mgry_u256(cbn_u256 const& n):
    n_(n)
  { }

  static mgry_u256 from_classical(cbn_u256 const& n) {
    return mgry_u256{mgry_mul(n, Rsq_p)};
  }

  mgry_u256 sqr() const {
    return mgry_u256{mgry_reduce(mul(n_,n_))};
  }

  cbn_u256 to_classical() const {
    const cbn_u512 wn = cbn::detail::pad<nlimbs>(n_);
    return mgry_reduce(wn);
  }

  cbn_u256 const& bn() const { return n_; }
  cbn_u256& bn() { return n_; }

  auto operator<=>(mgry_u256 const& o) const {
    return bn() <=> o.bn();
  }

private:
  cbn_u256 n_;
};

mgry_u256 operator*(mgry_u256 const& a, mgry_u256 const& b) {
  return mgry_u256{mgry_mul(a.bn(), b.bn())};
}

cbn_u256 rnd_u256() {
  cbn_u256 ret;
  getentropy(&ret, sizeof(ret));
  return ret;
}

} // namespace


int main()
{
  const cbn_u256 a = rnd_u256();
  const cbn_u256 b = rnd_u256();
  const auto ref = cbn::mul(a,b)%cbn::to_big_int(P);

  const auto ma = mgry_u256::from_classical(a);
  const auto mb = mgry_u256::from_classical(b);
  const auto mmul = ma*mb;
  const cbn_u256 mul = mmul.to_classical();

  std::cout << a << " + " << b << " % " << cbn::to_big_int(P) << std::endl;
  std::cout << ref << std::endl;
  std::cout << mul << std::endl;
  std::cout << (ref == mul) << std::endl;

  return 0;
}
