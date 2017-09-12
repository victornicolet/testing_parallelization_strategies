//
// Created by nicolet on 12/09/17.
//
#include <string>
#include <iostream>
using namespace std;

#ifndef PARALLEL_STRATEGIES_TESTING_UTILS_H
#define PARALLEL_STRATEGIES_TESTING_UTILS_H

struct result_data {
    double time_sequential;
    double time_strategy1;
    double time_strategy2;
    long pb_size;
    string name;
};

static void csvline(ostream& out, result_data rdat) {
    out << rdat.name << "," << rdat.pb_size << "," << rdat.time_strategy1 << "," << rdat.time_strategy2 << endl;
}

static double dsum (double *array, int length) {
    double _dsum = 0;
    for (int i = 0; i < length; i++)
        _dsum += array[i];
    return _dsum;
}

static double dmean (double* array, int length) {
    return  dsum(array, length) / length;
}

class Utils {

};


#endif //PARALLEL_STRATEGIES_TESTING_UTILS_H
