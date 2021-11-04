#include <ctbignum/bigint.hpp>
#include <ctbignum/mult.hpp>
#include <ctbignum/decimal_literals.hpp>
#include <ctbignum/montgomery.hpp>
#include <benchmark/benchmark.h>

#include <cstdint>
#include <array>
#include <random>
#include <iostream>

using namespace cbn::literals;

namespace {

std::random_device g_rnd;

template <class Bignum, bool LastZero = false>
auto random_bn() {
  using value_type = typename Bignum::value_type;
  Bignum ret;
  std::generate(std::begin(ret), std::end(ret), [&]() { return std::uniform_int_distribution<value_type>{}(g_rnd); });
  if constexpr (LastZero) {
    ret.back() &= 0xFF;
  }
  return ret;
}

template <class Bignum>
void bench_mul(benchmark::State& S) {
  const auto bn0 = random_bn<Bignum>();
  const auto bn1 = random_bn<Bignum>();

  for (auto _: S) {
    for (size_t i = 0; i < 4; ++i) {
      benchmark::DoNotOptimize(cbn::mul(bn0, bn1));
    }
  }
}

template <class Bignum>
void bench_add(benchmark::State& S) {
  const auto bn0 = random_bn<Bignum>();
  const auto bn1 = random_bn<Bignum>();

  for (auto _: S) {
    for (size_t i = 0; i < 4; ++i) {
      benchmark::DoNotOptimize(cbn::add(bn0, bn1));
    }
  }
}

template <class Bignum>
void bench_mgry_mul(benchmark::State& S) {
  const auto bn0 = random_bn<Bignum, true>();
  const auto bn1 = random_bn<Bignum, true>();

  for (auto _: S) {
    for (size_t i = 0; i < 4; ++i) {
      benchmark::DoNotOptimize(cbn::montgomery_mul(bn0, bn1,
            115792089237316195423570985008687907853269984665640564039457584007908834671663_Z));
    }
  }
}

} // anonymous

int main(int argc, char** argv)
{
  benchmark::RegisterBenchmark("add_256_u64_x4", &bench_add<cbn::big_int<4, uint64_t>>);

  benchmark::RegisterBenchmark("mul_128_u64_x4", &bench_mul<cbn::big_int<2, uint64_t>>);
  benchmark::RegisterBenchmark("mul_128_u32_x4", &bench_mul<cbn::big_int<4, uint32_t>>);

  benchmark::RegisterBenchmark("mul_256_u64_x4", &bench_mul<cbn::big_int<4, uint64_t>>);
  benchmark::RegisterBenchmark("mul_256_u32_x4", &bench_mul<cbn::big_int<8, uint32_t>>);

  benchmark::RegisterBenchmark("mgry_mul_256_u64_x4", &bench_mgry_mul<cbn::big_int<4, uint64_t>>);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
