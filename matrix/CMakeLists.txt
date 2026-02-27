cmake_minimum_required(VERSION 3.31)
project(matrix)

set(CMAKE_CXX_STANDARD 20)

add_executable(matrix main.cpp DenseMatrix.cpp DenseMatrix.hpp CSRMatrix.hpp CSRMatrix.hpp VectorOperations.cpp VectorOperations.hpp)

add_compile_options(-Wall -Wextra -Wpedantic -O2)