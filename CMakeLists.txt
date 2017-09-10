cmake_minimum_required(VERSION 2.8)

project(factormod)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -D_GLIBCXX_DEBUG")
set(main_SRC main.cpp poly.cpp matrix.cpp)

add_executable(main ${main_SRC})
