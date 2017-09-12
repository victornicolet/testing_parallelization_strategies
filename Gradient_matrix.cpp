//
// Created by nicolet on 12/09/17.
//

#include "Gradient_matrix.h"
#include "Stopwatch.h"
#include <tbb/tbb.h>

using namespace tbb;

struct GradMatrixPar {
    data_type **M;
    data_type *lower_border;
    data_type *top_border;
    iter_type b,e;
    bool ord;
    iter_type n;

    GradMatrixPar(data_type** _input, iter_type n) : ord(true), b(0), e(0), n(n), M(_input) {
        lower_border = new data_type[n];
        top_border = new data_type[n];
    }

    GradMatrixPar( GradMatrixPar& s, split ) {
        ord = true; b = 0; e = 0; n = s.n; M = s.M;
        lower_border = new data_type[n];
        top_border = new data_type[n];
    }

    void operator()( const blocked_range<iter_type>& r ) {
        data_type **a = M;
        for(iter_type i = r.begin(); i < r.end(); i++){
            for(iter_type j = 0; j < n; j++) {
                ord = ord 	&& a[i][j] > a[i+1][j]
                      && a[i][j] > a[i][j+1];
            }
        }

        b = r.begin();
        e = r.end() -1;

        for(iter_type j = 0; j < n; j++) {
            lower_border[j] = a[b][j];
            top_border[j] = a[e][j];
        }

    }
    void join( GradMatrixPar& rhs )
    {
        ord = ord && rhs.ord;
        for(iter_type i = 0; i < n; i++) {
            ord = ord && lower_border[i] > rhs.top_border[i];
        }
        delete rhs.top_border;
        delete lower_border;
        lower_border = rhs.lower_border;
    }

};

result_data Gradient_matrix::testGradientMatrix(data_type **input_matrix, iter_type pb_size) {
    StopWatch* t = new StopWatch();
    result_data pb_def = {0.0, 0.0, 0.0, pb_size, "Gradient_matrix"};
    double seq_time = 0.0;
    double tiled_version = 0.0;
    double split_loops_version = 0.0;

    // Sequential version
    bool ord = true;
    data_type ** a = input_matrix;
    t->start();
    for(iter_type i = 0; i < pb_size; i++){
        for(iter_type j = 0; j < pb_size; j++) {
            ord = ord 	&& a[i][j] > a[i+1][j]
                  && a[i][j] > a[i+1][j+1]
                  && a[i][j] > a[i][j+1];
        }
    }
    seq_time = t->stop();

    // Strategy 1 : with join updating using borders.
    GradMatrixPar gm(input_matrix, pb_size);
    t->clear();
    t->start();
    parallel_reduce(blocked_range<iter_type>(0,pb_size - 1), gm);
    tiled_version = t->stop();

    pb_def.time_sequential = seq_time;
    pb_def.time_strategy1 = tiled_version;
    pb_def.time_strategy2 = split_loops_version;

    cout << "Parallelization strategy 1 (horizontal stripes with border) : " << tiled_version << endl;
    cout << "Parallelization strategy 2 (split loops) : " <<  split_loops_version << endl;
    cout << "Sequential time : " << seq_time << endl;

    return pb_def;
}