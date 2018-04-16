#ifndef PARALLEL_STRATEGIES_TESTING_GRADIENT_MATRIX_H
#define PARALLEL_STRATEGIES_TESTING_GRADIENT_MATRIX_H


#include "../Stopwatch.h"
#include "../Utils.h"
#include <tbb/tbb.h>

class Gradient_matrix {
public:
    static result_data testGradientMatrix(data_type **, iter_type, test_params);
};


#endif //PARALLEL_STRATEGIES_TESTING_GRADIENT_MATRIX_H
