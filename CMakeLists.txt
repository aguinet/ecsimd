cmake_minimum_required(VERSION 3.21)
project(ecsimd)

enable_testing()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
	URL_HASH SHA256=353571c2440176ded91c2de6d6cd88ddd41401d14692ec1f99e35d013feda55a
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# TODO: bring eve using FetchContent

add_library(ecsimd INTERFACE)
target_include_directories(ecsimd INTERFACE ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/third-party)
set_target_properties(ecsimd PROPERTIES
	CXX_STANDARD 20
)

add_subdirectory(lib)
add_subdirectory(tests)
add_subdirectory(benchs)

add_executable(ecp256 src/ecp256.cpp)
target_include_directories(ecp256 PRIVATE ${PROJECT_SOURCE_DIR}/third-party)
