cmake_minimum_required(VERSION 3.5)
project(parallel_strategies_testing)

link_libraries(tbb)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ltbb")

set(SOURCE_FILES main.cpp Stopwatch.cpp Utils.cpp Utils.h Gradient_matrix.cpp Gradient_matrix.h Gradient_matrix_variations.cpp Gradient_matrix_variations.h ExamplesTaskBased.cpp ExamplesTaskBased.h)

add_executable(parallel_strategies_testing ${SOURCE_FILES})