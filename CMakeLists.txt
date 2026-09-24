cmake_minimum_required(VERSION 3.24)
project(CAPS VERSION 0.1.0 LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(CAPS_BUILD_TESTS "Build unit tests" ON)

# Compiler settings for every target: sources are UTF-8 (Å, ⟨⟩ in messages); MSVC needs to be told, and M_PI needs
# _USE_MATH_DEFINES there.
if(MSVC)
  add_compile_options(/utf-8 /W3 /bigobj /permissive-)
  add_compile_definitions(_USE_MATH_DEFINES NOMINMAX _CRT_SECURE_NO_WARNINGS)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()

# Fortran 2018 numeric kernels (k-space electrostatics, …): standard-conforming, optimised; the runtime is linked
# statically into the shipped library where the compiler allows it.
set(CMAKE_Fortran_STANDARD 2018)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2018 -Wall -fimplicit-none")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  # the kernels need no Fortran runtime (core/src/fortran_runtime.cpp provides the two entry points they can reach):
  # nothing of libgfortran or libquadmath is linked, so no platform has to ship them
  set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "")
  set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "")
endif()

# Windows builds use MinGW-w64 (gcc, g++, gfortran): link its C, C++ and thread runtimes statically so caps.dll and
# caps.exe need nothing beside them
if(MINGW)
  add_link_options(-static -static-libgcc -static-libstdc++)
  add_compile_definitions(_USE_MATH_DEFINES)   # M_PI under strict ISO C++
endif()

find_package(Threads REQUIRED)
# zlib: the system library (macOS, Linux); built from source where there is none (Windows)
find_package(ZLIB QUIET)
if(NOT ZLIB_FOUND)
  include(FetchContent)
  FetchContent_Declare(zlib
    URL https://github.com/madler/zlib/releases/download/v1.3.1/zlib-1.3.1.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  set(ZLIB_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(CMAKE_POLICY_VERSION_MINIMUM 3.5)   # zlib's CMakeLists still declares an older minimum
  FetchContent_MakeAvailable(zlib)
  add_library(ZLIB::ZLIB ALIAS zlibstatic)
  target_include_directories(zlibstatic INTERFACE ${zlib_SOURCE_DIR} ${zlib_BINARY_DIR})
endif()

add_subdirectory(core)
add_subdirectory(capi)
add_subdirectory(cli)

if(CAPS_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
