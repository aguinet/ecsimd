#include <ecsimd/bignum.h>
#include <ecsimd/add.h>
#include <ecsimd/mul.h>
#include <ecsimd/serialization.h>

#include <benchmark/benchmark.h>

#include <cstdint>
#include <random>

namespace {

std::random_device g_rnd;

template <class Bignum>
Bignum random_bn() {
  std::array<uint8_t, sizeof(Bignum)> ar;
  std::generate(std::begin(ar), std::end(ar), [&]() { return g_rnd(); });
  return ecsimd::bn_from_bytes_BE<Bignum>(ar);
}

template <class Bignum>
void bench_mul(benchmark::State& S) {
  ecsimd::wide_bignum<Bignum> bn0([](auto i, auto _) { return random_bn<Bignum>(); });
  ecsimd::wide_bignum<Bignum> bn1([](auto i, auto _) { return random_bn<Bignum>(); });

  for (auto _: S) {
    benchmark::DoNotOptimize(ecsimd::mul2(bn0, bn1));
  }
}

} // anonymous

int main(int argc, char** argv)
{
  benchmark::RegisterBenchmark("mul_128", &bench_mul<ecsimd::bignum_128>);
  benchmark::RegisterBenchmark("mul_256", &bench_mul<ecsimd::bignum_256>);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
