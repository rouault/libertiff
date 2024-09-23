cmake_minimum_required(VERSION 3.16)
project(libertiff_tests LANGUAGES CXX)
include(CTest)
set(CMAKE_CXX_STANDARD 11)

add_executable(demo demo.cpp)

add_subdirectory(tests)
