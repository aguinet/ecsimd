find_package(benchmark REQUIRED)

add_executable(bench_ops ops.cpp)
target_link_libraries(bench_ops ecsimd benchmark::benchmark)

add_executable(bench_ops_cbn ops_cbn.cpp)
target_link_libraries(bench_ops_cbn ecsimd benchmark::benchmark)

add_executable(bench_curve_group curve_group.cpp)
target_link_libraries(bench_curve_group ecsimd benchmark::benchmark)

add_executable(bench_p256_ref p256_ref.cpp)
target_include_directories(bench_p256_ref PRIVATE /usr/include/botan-2 /usr/include/crypto++)
target_link_libraries(bench_p256_ref benchmark::benchmark botan-2 cryptopp crypto)
