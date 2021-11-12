find_package(benchmark REQUIRED)

add_executable(bench_ops ops.cpp)
target_link_libraries(bench_ops ecsimd benchmark::benchmark)

add_executable(bench_ops_cbn ops_cbn.cpp)
target_link_libraries(bench_ops_cbn ecsimd benchmark::benchmark)

add_executable(bench_curve_group curve_group.cpp)
target_link_libraries(bench_curve_group ecsimd benchmark::benchmark)
