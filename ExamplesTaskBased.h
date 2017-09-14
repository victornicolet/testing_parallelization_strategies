//
// Created by nicolet on 13/09/17.
//

#ifndef PARALLEL_STRATEGIES_TESTING_EXAMPLESTASKBASED_H
#define PARALLEL_STRATEGIES_TESTING_EXAMPLESTASKBASED_H

#include "Utils.h"
#include "tbb/flow_graph.h"

result_data testMaxTopLeftSquareReduction(data_type**, iter_type, test_params);
result_data testMaxTopLeftSquareTaskPipelined(data_type**, iter_type, test_params);
#endif //PARALLEL_STRATEGIES_TESTING_EXAMPLESTASKBASED_H
