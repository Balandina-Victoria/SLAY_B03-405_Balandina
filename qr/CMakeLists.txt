cmake_minimum_required(VERSION 3.31)
project(qr)

set(CMAKE_CXX_STANDARD 20)

add_executable(qr DenseMatrix.cpp DenseMatrix.hpp hausholder.cpp)

add_compile_options(-Wall -Wextra -Wpedantic -O2)