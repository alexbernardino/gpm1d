cmake_minimum_required(VERSION 3.6)
project(gpm1d)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(gpm1d ${SOURCE_FILES})