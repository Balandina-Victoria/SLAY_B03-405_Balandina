cmake_minimum_required(VERSION 3.31)
project(solvers)

set(CMAKE_CXX_STANDARD 20)

add_executable(solvers main.cpp CSRMatrix.cpp VectorOperations.cpp)

add_compile_options(-Wall -Wextra -Wpedantic -O2)