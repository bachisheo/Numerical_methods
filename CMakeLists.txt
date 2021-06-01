cmake_minimum_required(VERSION 3.19)
project(Numerical_methods)

set(CMAKE_CXX_STANDARD 14)

add_executable(Numerical_methods main.cpp num_methods.cpp num_methods.h vect.cpp vect.h matrix.cpp matrix.h geometry.cpp geometry.h)