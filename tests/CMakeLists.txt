include(GoogleTest)

add_executable(ops ops.cpp)
target_link_libraries(ops ecsimd gtest_main)

add_executable(mgry mgry.cpp)
target_link_libraries(mgry ecsimd gtest_main)

add_executable(curve_point curve_point.cpp)
target_link_libraries(curve_point ecsimd gtest_main)

add_executable(curve_group curve_group.cpp)
target_link_libraries(curve_group ecsimd gtest_main)

gtest_discover_tests(ops)
gtest_discover_tests(mgry)
gtest_discover_tests(curve_point)
gtest_discover_tests(curve_group)
