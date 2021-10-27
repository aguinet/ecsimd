enable_testing()

add_executable(ops ops.cpp)
target_link_libraries(ops ecsimd gtest_main)

include(GoogleTest)
gtest_discover_tests(ops)
