cmake_minimum_required(VERSION 3.5)
project(parallel_strategies_testing)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set(SOURCE_FILES MaxBlockOnesCustomSplit.cpp)

add_executable(mbo_test ${SOURCE_FILES})