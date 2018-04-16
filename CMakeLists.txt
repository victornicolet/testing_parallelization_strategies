cmake_minimum_required(VERSION 3.5)
project(parallel_strategies_testing)

link_libraries(tbb)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ltbb -fopenmp")

set(SOURCE_FILES main.cpp Stopwatch.cpp Utils.cpp Utils.h Gradient_matrix.cpp
        Gradient_matrix.h Gradient_matrix_variations.cpp Gradient_matrix_variations.h
        ExamplesTaskBased.cpp ExamplesTaskBased.h Image.cpp Image.h MaxTopStrip.cpp
        MaxTopStrip.h TestDepFlowGraph.cpp TestDepFlowGraph.h
        longest_common_subsequence.cpp longest_common_subsequence.h MaxTopRightRectangle.cpp MaxTopRightRectangle.h)

OPTION(DEFINE_MTRR_LWD_AUX "Compile wrong version with leftwards auxiliary for MTRR" OFF)

IF(DEFINE_MTRR_LWD_AUX)
    add_definitions(-DMTRR_LEFTWARDS_AUX)
ENDIF(DEFINE_MTRR_LWD_AUX)

add_executable(parallel_strategies_testing ${SOURCE_FILES})