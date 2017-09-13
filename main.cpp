#define UT_THREAD_DEFAULT_STACK_SIZE (8U*1024U*1024U)


#include <iostream>
#include <unistd.h>
#include <fstream>
#include "Stopwatch.h"
#include "Utils.h"
#include "Gradient_matrix.h"
#include "Gradient_matrix_variations.h"
#include "ExamplesTaskBased.h"
#include <tbb/tbb.h>

#define TESTS_NUM 40
#define SIZES 10

using namespace tbb;
using namespace std;

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
result_data topMaxStripStrategyComparison(data_type** input, result_data pb_def ) {
    StopWatch *t = new StopWatch;
    iter_type n = pb_def.pb_size;
    // Try the row-chunk strategy first (outer loop in parallel)
    double* row_par_time = new double[TESTS_NUM];
    if (input) {
        MaxStripPar mstrip_r(n, input);
        for(int i = 0; i < TESTS_NUM; i++) {
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
    double* inner_strategy = new double[TESTS_NUM];
    t->clear();
    data_type mtops = 0;
    data_type strip_sum = 0;
    data_type topsum = 0;
    SumRow sumrow(n);
    for(int test = 0; test < TESTS_NUM; test++) {
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

    double row_par_time_mean = dmean(row_par_time, TESTS_NUM);
    double inner_strategy_mean = dmean(inner_strategy, TESTS_NUM);
    cout << "Parallelization strategy 1 (horizontal stripes) : " << row_par_time_mean << endl;
    cout << "Parallelization strategy 2 (inner loop) : " << inner_strategy_mean << endl;
    cout << "Sequential time : " << time_seq << endl;
    pb_def.time_strategy1 = row_par_time_mean;
    pb_def.time_strategy2 = inner_strategy_mean;
    pb_def.time_sequential = time_seq;
    return pb_def;
}

result_data testMaxStripStrategyComparison(iter_type pb_size) {
    result_data pb_def = { 0.0, 0.0, 0.0, pb_size, "max_top_strip"};

    data_type** _data = init_data_matrix(pb_size);
    cout << "Test different strategies for maxtopstrip ..." << endl;
    result_data pb_res = topMaxStripStrategyComparison(_data, pb_def);
    clean_data_matrix(_data, pb_size);

    return pb_res;
}

result_data testGradientMatrixStrategyComparison(iter_type pb_size,
                                                 data_type** _data,
                                                 test_params tp) {
    cout << "Test different strategies for gradientmatrix  ..." << endl;
    return Gradient_matrix::testGradientMatrix(_data, pb_size, tp);
}


test_data testGradientMatrixStrategyComparison2(iter_type pb_size, data_type** _data,
test_params tp) {
    cout << "Test different strategies for gradientmatrix (var) ..." << endl;
    result_data pb_res1 =
            Gradient_matrix_variations::testGradientMatrixVariation1(_data, pb_size, tp);
    result_data pb_res2 =
            Gradient_matrix_variations::testGradientMatrixVariation2(_data, pb_size, tp);
    result_data* pb_res = new result_data[2];
    pb_res[0] = pb_res1;
    pb_res[1] = pb_res2;
    return {pb_res, 2};
}

void run_testset1(test_params tp) {
    tp.out << tp.test_names << endl;
    for (int i = 0; i < tp.nsizes; ++i) {
        data_type** m = init_data_matrix(tp.pb_sizes[i]);
        result_data r = testGradientMatrixStrategyComparison(tp.pb_sizes[i], m, tp);
        test_data pb_ress = testGradientMatrixStrategyComparison2(tp.pb_sizes[i],m,tp);
//      "tname,seq,single,split,var1,var2"
        tp.out << tp.pb_sizes[i] << ","
               << r.time_sequential << ","
               << r.time_strategy1 << ","
               << r.time_strategy2 << ","
               << pb_ress.results[0].time_strategy1 << ","
               << pb_ress.results[1].time_strategy1 << endl;
        clean_data_matrix(m, tp.pb_sizes[i]);
    }
}

void run_testset2(test_params tp) {
    tp.out << tp.test_names << endl;
    for (int i = 0; i < tp.nsizes; ++i) {
        data_type** m = init_data_matrix(tp.pb_sizes[i]);
        cout << "Speedup : " <<
             testMaxTopLeftSquareReduction(m,tp.pb_sizes[i],tp) << endl;
        clean_data_matrix(m, tp.pb_sizes[i]);
    }
}

int main(int argc, char** argv) {
    int num_cores = -1;
    int opt;

    while ((opt = getopt(argc, argv, "n:")) != -1) {
        switch (opt) {
            case 'n':
                num_cores = atoi(optarg);
                break;
            default: /* '?' */
                cerr << "Usage: " << argv[0] << " [-n ncores]\n";
                exit(EXIT_FAILURE);
        }
    }

    cout << "Number of cores : " << num_cores << endl;

    // If a number of cores is specified and greate than 1,
    // then set then number of cores used by TBB.
    if (num_cores > 1) {
        static task_scheduler_init
        init(task_scheduler_init::deferred);

        init.initialize(num_cores, UT_THREAD_DEFAULT_STACK_SIZE);
    }

    // Test for matrix sizes ...
    int nsizes = 20;
    iter_type* pb_sizes = new iter_type[nsizes];

    for (int j = 0; j < nsizes; ++j) {
        pb_sizes[j] = (1 << 15) + 10000 * j;
    }



    ofstream out_csv;
    out_csv.open("parallel_strategies.csv");


//    test_params test_params1 = {
//            nsizes,
//            pb_sizes,
//            "size,seq,single,split,var1,var2",
//            TESTS_NUM,
//            out_csv
//    };
//
//    run_testset1(test_params1);

    test_params test_params2 = {
            nsizes,
            pb_sizes,
            "speedup parallel_reduce/sequential",
            TESTS_NUM,
            out_csv
    };

    run_testset2(test_params2);

    out_csv.close();
    delete pb_sizes;
    return 0;
}