#include <ecsimd/jacobian_curve_point.h>
#include <ecsimd/curve_group.h>
#include <ecsimd/curve_nist_p256.h>
#include <ecsimd/serialization.h>
#include <ecsimd/literals.h>

#include <benchmark/benchmark.h>

#include <cstdint>
#include <random>

using namespace ecsimd;
using namespace ecsimd::literals;

namespace {

template <concepts::wide_bignum WBN>
auto wide_bignum_set1(std::array<uint8_t, sizeof(typename WBN::value_type)> const& bytes) {
  const auto BN = bn_from_bytes_BE<typename WBN::value_type>(bytes);
  return WBN{BN};
}

void bench_p256(benchmark::State& S) {
  using Curve = curve_nist_p256;
  using CurveGroup = curve_group<Curve>;
  using WBN = curve_wide_bn_t<Curve>;

  const auto WJG = CurveGroup::WJG();
  const auto x = wide_bignum_set1<WBN>("0a891cecc2bf13b0aca744434a9c9f4bd7bf5c8ed86e2f76e7df72bad813bd80"_hex);

  for (auto _: S) {
    const auto WMP = CurveGroup::scalar_mult(x, WJG);
    benchmark::DoNotOptimize(WMP.to_affine());
  }
}
void bench_p256_1s(benchmark::State& S) {
  using Curve = curve_nist_p256;
  using CurveGroup = curve_group<Curve>;
  using WBN = curve_wide_bn_t<Curve>;

  const auto WJG = CurveGroup::WJG();
  const auto x = bn_from_bytes_BE<typename WBN::value_type>("0a891cecc2bf13b0aca744434a9c9f4bd7bf5c8ed86e2f76e7df72bad813bd80"_hex);

  for (auto _: S) {
    const auto WMP = CurveGroup::scalar_mult_1s(x, WJG);
    benchmark::DoNotOptimize(WMP.to_affine());
  }
}

} // anonymous

int main(int argc, char** argv)
{
  benchmark::RegisterBenchmark("scalar_mult_p256_x4", bench_p256)->Unit(benchmark::kMicrosecond);
  benchmark::RegisterBenchmark("scalar_mult_p256_1s_x4", bench_p256_1s)->Unit(benchmark::kMicrosecond);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
