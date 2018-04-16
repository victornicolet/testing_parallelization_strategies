cmake_minimum_required(VERSION 3.5)
project(parallel_strategies_testing)

link_libraries(tbb)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ltbb -fopenmp")

set(SOURCE_FILES main.cpp Stopwatch.cpp Utils.cpp Utils.h n2/Gradient_matrix.cpp
        n2/Gradient_matrix.h n2/Gradient_matrix_variations.cpp n2/Gradient_matrix_variations.h
        ExamplesTaskBased.cpp ExamplesTaskBased.h Image.cpp Image.h n2/MaxTopStrip.cpp
        n2/MaxTopStrip.h TestDepFlowGraph.cpp TestDepFlowGraph.h
        n2/longest_common_subsequence.cpp n2/longest_common_subsequence.h n2/MaxTopRightRectangle.cpp n2/MaxTopRightRectangle.h)

OPTION(DEFINE_MTRR_LWD_AUX "Compile wrong version with leftwards auxiliary for MTRR" OFF)

IF(DEFINE_MTRR_LWD_AUX)
    add_definitions(-DMTRR_LEFTWARDS_AUX)
ENDIF(DEFINE_MTRR_LWD_AUX)

add_executable(parallel_strategies_testing ${SOURCE_FILES})