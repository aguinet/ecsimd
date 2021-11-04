#include <ecsimd/bignum.h>
#include <ecsimd/add.h>
#include <ecsimd/mul.h>
#include <ecsimd/mgry_mul.h>
#include <ecsimd/serialization.h>
#include <ecsimd/literals.h>

#include <benchmark/benchmark.h>

#include <cstdint>
#include <random>

using namespace ecsimd;
using namespace ecsimd::literals;

namespace {

std::random_device g_rnd;

struct P {
  static constexpr auto value = bn_from_bytes_BE<bignum_256>("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F"_hex);
};

template <class Bignum, bool LastZero = false>
Bignum random_bn() {
  std::array<uint8_t, sizeof(Bignum)> ar;
  std::generate(std::begin(ar), std::end(ar), [&]() { return g_rnd(); });
  if constexpr (LastZero) {
    ar[sizeof(Bignum)-1] = 0;
  }
  return bn_from_bytes_BE<Bignum>(ar);
}

template <class Bignum>
void bench_mul(benchmark::State& S) {
  wide_bignum<Bignum> bn0([](auto i, auto _) { return random_bn<Bignum>(); });
  wide_bignum<Bignum> bn1([](auto i, auto _) { return random_bn<Bignum>(); });

  for (auto _: S) {
    benchmark::DoNotOptimize(mul(bn0, bn1));
  }
}

template <class Bignum>
void bench_add(benchmark::State& S) {
  wide_bignum<Bignum> bn0([](auto i, auto _) { return random_bn<Bignum>(); });
  wide_bignum<Bignum> bn1([](auto i, auto _) { return random_bn<Bignum>(); });

  for (auto _: S) {
    benchmark::DoNotOptimize(add(bn0, bn1));
  }
}

void bench_mgry_mul(benchmark::State& S) {
  using BN = bignum_256;
  wide_bignum<BN> bn0([](auto i, auto _) { return random_bn<BN, true>(); });
  wide_bignum<BN> bn1([](auto i, auto _) { return random_bn<BN, true>(); });

  for (auto _: S) {
    benchmark::DoNotOptimize(details::mgry_mul<P>(bn0, bn1));
  }
}

} // anonymous

int main(int argc, char** argv)
{
  benchmark::RegisterBenchmark("add_256", &bench_add<bignum_128>);

  benchmark::RegisterBenchmark("mul_128", &bench_mul<bignum_128>);
  benchmark::RegisterBenchmark("mul_256", &bench_mul<bignum_256>);

  benchmark::RegisterBenchmark("mgry_mul_256", bench_mgry_mul);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}