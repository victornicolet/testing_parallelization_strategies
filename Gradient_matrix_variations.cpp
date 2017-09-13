//
// Created by nicolet on 13/09/17.
//

#include "Gradient_matrix_variations.h"
#include "Stopwatch.h"
#include <tbb/tbb.h>

using namespace tbb;

struct GradMatrixVar1Par {
    data_type **M;
    iter_type b,e;
    bool ord;
    iter_type n;

    GradMatrixVar1Par(data_type** _input, iter_type n) :
            ord(true), b(0), e(0), n(n), M(_input) {}

    GradMatrixVar1Par( GradMatrixVar1Par& s, split ) {
        ord = true; b = 0; e = 0; n = s.n; M = s.M;
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

    }
    void join( GradMatrixVar1Par& rhs )
    {
        ord = ord && rhs.ord;
        // Doing useless linear computation, for testing purposese only
        for(iter_type i = 0; i < n; i++) {
            ord = ord && (1 << i) > i;
        }
    }

};

struct GradMatrixVar1Row {
    data_type **M;
    bool ord;
    iter_type n;

    GradMatrixVar1Row(data_type** _input, iter_type n) : ord(true), M(_input), n(n) {}
    GradMatrixVar1Row(GradMatrixVar1Row& s, split): ord(true), M(s.M), n(s.n) {}

    void operator()(const blocked_range<iter_type>& r){
        data_type** a = M;
        for(iter_type i = r.begin(); i != r.end(); i++){
            for(iter_type j = 0; j < n - 1; j++) {
                ord = ord && a[i][j] < a[i][j+1];
            }
        }
    }

    void join(GradMatrixVar1Row& rhs) {
        ord = ord && rhs.ord;
    }
};

struct GradMatrixVar1Col {
    data_type **M;
    bool ord;
    iter_type n;

    GradMatrixVar1Col(data_type** _input, iter_type n) : ord(true), M(_input), n(n) {}
    GradMatrixVar1Col(GradMatrixVar1Col& s, split): ord(true), M(s.M), n(s.n) {}

    void operator()(const blocked_range<iter_type>& r){
        data_type** a = M;
        for(iter_type i = r.begin(); i != r.end(); i++){
            for(iter_type j = 0; j < n - 1; j++) {
                // Inverted row/col indices
                ord = ord && a[j][i] < a[j+1][i];
            }
        }
    }

    void join(GradMatrixVar1Col& rhs) {
        ord = ord && rhs.ord;
    }
};


result_data Gradient_matrix_variations::testGradientMatrixVariation1(
        data_type **input_matrix,
        iter_type pb_size) {
    StopWatch* t = new StopWatch();
    result_data pb_def = {0.0, 0.0, 0.0, pb_size, "GM_linearTimeJoin_constantMem"};
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
    GradMatrixVar1Par gm(input_matrix, pb_size);
    t->clear();
    t->start();
    parallel_reduce(blocked_range<iter_type>(0,pb_size - 1), gm);
    tiled_version = t->stop();


    // Strategy 2 : split the loop in two steps: check columns are orderd and lines
    GradMatrixVar1Row gr(input_matrix, pb_size);
    GradMatrixVar1Col gc(input_matrix, pb_size);
    t->clear();
    t->start();
    parallel_reduce(blocked_range<iter_type>(0,pb_size - 1), gr);
    parallel_reduce(blocked_range<iter_type>(0,pb_size - 1), gc);
    // Don't forget to join the two results
    bool ord_split = gr.ord && gc.ord;
    split_loops_version = t->stop();

    pb_def.time_sequential = seq_time;
    pb_def.time_strategy1 = tiled_version;
    pb_def.time_strategy2 = split_loops_version;

    cout << "Parallelization strategy 1 (horizontal stripes with border) : "
         << tiled_version << endl;
    cout << "Parallelization strategy 2 (split loops) : "
         <<  split_loops_version << endl;
    cout << "Sequential time : " << seq_time << endl;

    return pb_def;
}


// Second variation

struct GradMatrixVar2Par {
    data_type **M;
    data_type *lower_border;
    data_type *top_border;
    iter_type b,e;
    bool ord;
    iter_type n;

    GradMatrixVar2Par(data_type** _input, iter_type n) : ord(true), b(0), e(0), n(n), M(_input) {
        lower_border = new data_type[n];
        top_border = new data_type[n];
    }

    GradMatrixVar2Par( GradMatrixVar2Par& s, split ) {
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

    // In this version the join does not perform any "linear" computation
    // but we are still passing the top and lower borders around.

    void join( GradMatrixVar2Par& rhs )
    {
        ord = ord && rhs.ord;
        delete rhs.top_border;
        delete lower_border;
        lower_border = rhs.lower_border;
    }

};


result_data Gradient_matrix_variations::testGradientMatrixVariation2(
        data_type **input_matrix,
        iter_type pb_size) {
    StopWatch* t = new StopWatch();
    result_data pb_def = {0.0, 0.0, 0.0, pb_size, "GM_constantTimeJoin_linearMem"};
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
    GradMatrixVar2Par gm(input_matrix, pb_size);
    t->clear();
    t->start();
    parallel_reduce(blocked_range<iter_type>(0,pb_size - 1), gm);
    tiled_version = t->stop();


    // Strategy 2 : split the loop in two steps: check columns are orderd and lines
    GradMatrixVar1Row gr(input_matrix, pb_size);
    GradMatrixVar1Col gc(input_matrix, pb_size);
    t->clear();
    t->start();
    parallel_reduce(blocked_range<iter_type>(0,pb_size - 1), gr);
    parallel_reduce(blocked_range<iter_type>(0,pb_size - 1), gc);
    // Don't forget to join the two results
    bool ord_split = gr.ord && gc.ord;
    split_loops_version = t->stop();

    pb_def.time_sequential = seq_time;
    pb_def.time_strategy1 = tiled_version;
    pb_def.time_strategy2 = split_loops_version;

    cout << "Parallelization strategy 1 (horizontal stripes with border) : "
         << tiled_version << endl;
    cout << "Parallelization strategy 2 (split loops) : "
         <<  split_loops_version << endl;
    cout << "Sequential time : " << seq_time << endl;

    return pb_def;
}
