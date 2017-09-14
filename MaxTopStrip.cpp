//
// Created by nicolet on 14/09/17.
//

#include "MaxTopStrip.h"
#include "Utils.h"
#include <tbb/tbb.h>
#include "Stopwatch.h"

// Max top strip. Parallelization strategy : Chunks of several rows
struct MaxStripPar {
    data_type **input;

    data_type mtops;
    data_type topsum;
    iter_type row_len;

    MaxStripPar(iter_type rl, data_type** _input) : mtops(0), topsum(0), row_len(rl), input(_input){}

    MaxStripPar( MaxStripPar& s, split ) {mtops = 0; topsum = 0; row_len = s.row_len; input = s.input;}

    void operator()( const blocked_range<iter_type>& r ) {
        data_type stripsum = 0;
        for(iter_type i = r.begin(); i != r.end(); ++i) {
            stripsum = 0;
            for(iter_type j = 0; j < row_len; j++) {
                stripsum += input[i][j];
            }
            topsum += stripsum;
            mtops = std::max(topsum, mtops);
        }
    }
    void join( MaxStripPar& rhs ) {mtops += std::max(mtops, topsum + rhs.mtops);}

};


// Max trop strip, inner loop (simply a sum tbh)
struct SumRow {
    data_type *input_row;
    data_type sum;
    iter_type row_len;

    SumRow(iter_type rl):sum(0), row_len(rl){}

    SumRow(SumRow& s, split) {sum = 0; row_len = s.row_len; input_row = s.input_row; }

    void operator()(const blocked_range<iter_type>& r) {
        for(iter_type i = r.begin(); i!= r.end(); ++i) {
            sum += input_row[i];
        }
    }

    void join(SumRow& rhs) {sum += rhs.sum; }
};



// Note: Reads A[0..n] and writes output[1..n-1].
result_data topMaxStripStrategyComparison(data_type** input, result_data pb_def, test_params tp) {
    StopWatch *t = new StopWatch;
    iter_type n = pb_def.pb_size;
    // Try the row-chunk strategy first (outer loop in parallel)
    double* row_par_time = new double[tp.number_per_test];
    if (input) {
        MaxStripPar mstrip_r(n, input);
        for(int i = 0; i < tp.number_per_test; i++) {
            t->start();
            parallel_reduce(blocked_range<iter_type>(1, n), mstrip_r);
            row_par_time[i] = t->stop();
            t->clear();
        }
    } else {
        cout << "No data." << endl;
        return pb_def;
    }
    // Inner loop in parallel strategy
    double* inner_strategy = new double[tp.number_per_test];
    t->clear();
    data_type mtops = 0;
    data_type strip_sum = 0;
    data_type topsum = 0;
    SumRow sumrow(n);
    for(int test = 0; test < tp.number_per_test; test++) {
        t->start();
        for (iter_type i = 0; i < n; i++) {
            sumrow.input_row = input[i];
            parallel_reduce(blocked_range<iter_type>(1, n), sumrow);
            topsum += sumrow.sum;
            mtops = std::max(mtops, topsum);
        }
        inner_strategy[test] = t->stop();
    }


    // Sequential time, for reference
    mtops = 0;
    topsum = 0;

    t->clear();
    t->start();

    for (iter_type i = 0; i < n; i++) {
        int sumr = 0;
        for(iter_type j = 0; j < n; j++) {
            sumr += input[i][j];
        }
        topsum += sumr;
        mtops = std::max(mtops, topsum);
    }
    double time_seq = t->stop();

    double row_par_time_mean = dmean(row_par_time, tp.number_per_test);
    double inner_strategy_mean = dmean(inner_strategy, tp.number_per_test);
    cout << "Parallelization strategy 1 (horizontal stripes) : " << row_par_time_mean << endl;
    cout << "Parallelization strategy 2 (inner loop) : " << inner_strategy_mean << endl;
    cout << "Sequential time : " << time_seq << endl;
    pb_def.time_strategy1 = row_par_time_mean;
    pb_def.time_strategy2 = inner_strategy_mean;
    pb_def.time_sequential = time_seq;
    return pb_def;
}