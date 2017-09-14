#ifndef PARALLEL_STRATEGIES_TESTING_UTILS_H
#define PARALLEL_STRATEGIES_TESTING_UTILS_H

#include <string>
#include <iostream>
#include <tbb/tbb.h>
#include "Image.h"

using namespace tbb;
using namespace std;

typedef long iter_type;
typedef int data_type;

typedef struct {
    int nsizes;
    iter_type* pb_sizes;
    string test_names;
    int number_per_test;
    ostream& out;
} test_params;

struct result_data {
    double time_sequential;
    double time_strategy1;
    double time_strategy2;
    long pb_size;
    string name;
};

typedef struct {result_data* results; int num_results;} test_data;

static void csvline(ostream& out, result_data rdat) {
    out << rdat.name << "," << rdat.pb_size << "," << rdat.time_sequential << "," << rdat.time_strategy1 << "," << rdat.time_strategy2 << endl;
}

static void csvlines(ostream& out, test_data t) {
    for (int i = 0; i < t.num_results; ++i) {
        csvline(out, t.results[i]);
    }
}

static double dsum (double *array, int length) {
    double _dsum = 0;
    for (int i = 0; i < length; i++)
        _dsum += array[i];
    return _dsum;
}

struct Initializer {
    data_type **M;
    iter_type pb_size;

    Initializer(data_type **_m, iter_type _s) : M(_m) , pb_size(_s){}

    void operator()(const blocked_range<iter_type>& r) const {
        for(iter_type i = r.begin(); i < r.end(); i++) {
            M[i] = new data_type[pb_size];
            for(iter_type j = 0; j < pb_size; j++) {
                M[i][j] = (((data_type) rand()) % 200) - 100;
            }
        }
    }
};

static data_type** init_data_matrix(iter_type pb_size) {
    cout << "Initialize data (size " << pb_size << ") ..." << endl;
    data_type** _data = new data_type*[pb_size];
    Initializer init(_data, pb_size);
    parallel_for(blocked_range<iter_type>(0,pb_size), init);
    return _data;
}

static data_type** init_data_matrix_seq(iter_type pb_size) {
    cout << "Initialize data (size " << pb_size << ") ..." << endl;
    data_type** M = new data_type*[pb_size];
    for(iter_type i = 0; i < pb_size; i++) {
        M[i] = new data_type[pb_size];
        for(iter_type j = 0; j < pb_size; j++) {
            M[i][j] = (((data_type) rand()) % 200) - 100;
        }
    }
    return M;
}

static void init_data_from_image(const char* filename){
    MImage img(filename);
    unchar** data = img.getData();
}

static void clean_data_matrix(data_type **_data, iter_type pb_size) {
    for(iter_type i = 0; i < pb_size; i++) {
        delete _data[i];
    }
    delete _data;
}

static data_type* mzeroes(iter_type n) {
    data_type *z = new data_type[n];
    for (iter_type i = 0; i < n; ++i) {
        z[i] = 0;
    }
    return z;
}

static double dmean (double* array, int length) {
    return  dsum(array, length) / length;
}

class Utils {

};


#endif //PARALLEL_STRATEGIES_TESTING_UTILS_H
