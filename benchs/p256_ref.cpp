// Botan
#include <botan/system_rng.h>
#include <botan/nist_keywrap.h>
#include <botan/ec_group.h>
#include <botan/pubkey.h>

// Crypto++
#include <integer.h>
#include <eccrypto.h>
#include <osrng.h>
#include <oids.h>

// OpenSSL
#include <openssl/obj_mac.h>
#include <openssl/ec.h>
#include <openssl/bn.h>

#include <benchmark/benchmark.h>

static void bench_botan(benchmark::State& St) {
  std::unique_ptr<Botan::RandomNumberGenerator> rng(new Botan::System_RNG{});

  Botan::EC_Group nist256p("secp256r1");
  Botan::PointGFp const& BP = nist256p.get_base_point();

  Botan::BigInt const x = nist256p.random_scalar(*rng);
  Botan::PointGFp p;

  for (auto _: St) {
    p = BP;
    p *= x;
    benchmark::DoNotOptimize(p);
  }
}

static void bench_cryptopp(benchmark::State& St) {
  using namespace CryptoPP;
  typedef DL_GroupParameters_EC<ECP> GroupParameters;
  typedef DL_GroupParameters_EC<ECP>::Element Element;

  AutoSeededRandomPool prng;
  GroupParameters group;
  group.Initialize(ASN1::secp256r1());

  // private key
  Integer const x(prng, Integer::One(), group.GetMaxExponent());

  for (auto _: St) {
    auto p = group.ExponentiateBase(x);
    benchmark::DoNotOptimize(p);
  }
}

static void bench_openssl(benchmark::State& St) {
  EC_GROUP *curve = EC_GROUP_new_by_curve_name(NID_X9_62_prime256v1);
  if (!curve) {
    printf("error: unable to load openssl curve!\n");
    return;
  }

  BN_CTX *ctx = BN_CTX_new();
  if (!ctx) {
    printf("error: unable to create an openssl BN ctx!\n");
    return;
  }
  BIGNUM *prv = BN_new();
  BIGNUM *order = BN_new();
  if (!prv || !order) {
    printf("error: unable to create an openssl BN object!\n");
    return;
  }
  if (!EC_GROUP_get_order(curve, order, ctx)) {
    printf("error: EC_GROUP_get_order\n");
    return;
  }
  if (!BN_rand_range(prv, order)) {
    printf("error: BN_rand_range_ex\n");
    return;
  }

  EC_POINT* pub = EC_POINT_new(curve);
  if (!EC_POINT_mul(curve, pub, prv, NULL, NULL, ctx)) {
    printf("error: EC_POINT_mul\n");
    return;
  }

  for (auto _: St) {
    EC_POINT_mul(curve, pub, prv, NULL, NULL, ctx);
  }
  EC_POINT_free(pub);
  BN_free(prv);
  BN_free(order);
  BN_CTX_free(ctx);
}

BENCHMARK(bench_botan)->Unit(benchmark::kMicrosecond);
BENCHMARK(bench_cryptopp)->Unit(benchmark::kMicrosecond);
BENCHMARK(bench_openssl)->Unit(benchmark::kMicrosecond);
BENCHMARK_MAIN();
