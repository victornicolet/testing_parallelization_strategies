//
// Created by nicolet on 12/09/17.
//
#include <string>
#include <iostream>
using namespace std;

#ifndef PARALLEL_STRATEGIES_TESTING_UTILS_H
#define PARALLEL_STRATEGIES_TESTING_UTILS_H

typedef long iter_type;
typedef int data_type;


struct result_data {
    double time_sequential;
    double time_strategy1;
    double time_strategy2;
    long pb_size;
    string name;
};

static void csvline(ostream& out, result_data rdat) {
    out << rdat.name << "," << rdat.pb_size << "," << rdat.time_sequential << "," << rdat.time_strategy1 << "," << rdat.time_strategy2 << endl;
}

static double dsum (double *array, int length) {
    double _dsum = 0;
    for (int i = 0; i < length; i++)
        _dsum += array[i];
    return _dsum;
}

static void init_data_matrix(data_type** _data, iter_type pb_size) {
    cout << "Initialize data (size " << pb_size << ") ..." << endl;
    for(iter_type i = 0; i < pb_size; i++) {
        _data[i] = new data_type[pb_size];
        for(iter_type j = 0; j < pb_size; j++) {
            _data[i][j] = (((data_type) rand()) % 200) - 100;
        }
    }
}

static void clean_data_matrix(data_type **_data, iter_type pb_size) {
    for(iter_type i = 0; i < pb_size; i++) {
        delete _data[i];
    }
    delete _data;
}

static double dmean (double* array, int length) {
    return  dsum(array, length) / length;
}

class Utils {

};


#endif //PARALLEL_STRATEGIES_TESTING_UTILS_H
