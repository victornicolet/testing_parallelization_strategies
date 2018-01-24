#define UT_THREAD_DEFAULT_STACK_SIZE (8U*1024U*1024U)


#include <iostream>
#include <unistd.h>
#include <fstream>
#include "Utils.h"
#include "MaxTopStrip.h"
#include "Gradient_matrix.h"
#include "Gradient_matrix_variations.h"
#include "ExamplesTaskBased.h"
#include "longest_common_subsequence.h"
#include "TestDepFlowGraph.h"
#define TESTS_NUM 40


using namespace tbb;
using namespace std;

result_data testMaxStripStrategyComparison(iter_type pb_size, test_params tp) {
    result_data pb_def = { 0.0, 0.0, 0.0, pb_size, "max_top_strip"};
    data_type** _data = init_data_matrix_seq(pb_size);
    cout << "Test different strategies for maxtopstrip ..." << endl;
    result_data pb_res = topMaxStripStrategyComparison(_data, pb_def, tp);
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

// Test set 1 : comparing strategies for sorted-matrix (gradient) (with sorted-array)
void run_testset1(bool sorted, test_params tp) {
    tp.out << tp.test_names << endl;
    for (int i = 0; i < tp.nsizes; ++i) {
        data_type** m = sorted ? init_data_matrix_sorted(tp.pb_sizes[i]) : init_data_matrix(tp.pb_sizes[i]);
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

// Test set 2 : comparing strategies for  max top left square
void run_testset2(test_params tp) {
    tp.out << tp.test_names << endl;
    for (int i = 0; i < tp.nsizes; ++i) {
        data_type** m = init_data_matrix_seq(tp.pb_sizes[i]);
        csvline(tp.out, testMaxTopLeftSquareReduction(m,tp.pb_sizes[i],tp));
//       csvline(tp.out, testMaxTopLeftSquareTaskPipelined(m,tp.pb_sizes[i],tp));
        csvline(tp.out, testMtlsMultiscan(m, tp.pb_sizes[i], tp));
        clean_data_matrix(m, tp.pb_sizes[i]);
    }
}

// Test set 3 : measuring performance of tiled parallel implementation of LCS
void run_testset3(test_params tp){
    tp.out << tp.test_names << endl;
    for (int i = 0; i < tp.nsizes; ++i) {
        long pbsize = tp.pb_sizes[i];
        LCS lcs(pbsize, pbsize);
        lcs.measure_perfs(10);
        LCS_timedata tmd = lcs.get_perfs();
        csvline(tp.out, tmd);
    }
}

int _main(int argc, char** argv) {
    int num_cores = -1;
    int start_pow2_size = 15;
    int opt;
    bool test_lcs = false;
    bool test_1 = false;
    bool test_2 = false;

    while ((opt = getopt(argc, argv, "n:s:lm")) != -1) {
        switch (opt) {
            case 'n':
                num_cores = atoi(optarg);
                break;
            case 's':
                start_pow2_size = atoi(optarg);
                break;
            case 'l':
                test_lcs = true;
                break;
            case 'm':
                test_1 = true;
                test_2 = true;
                break;
            default: /* '?' */
                cerr << "Usage: " << argv[0] << " [-n ncores] [-s power of 2 size]\n";
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
        pb_sizes[j] = (1 << start_pow2_size) + 10000 * j;
    }



    if(test_1) {
        ofstream out_csv;
        out_csv.open("parallel_strategies_test_sortedness.csv");


        test_params test_params1 = {
                nsizes,
                pb_sizes,
                "size,seq,single,split,var1,var2",
                TESTS_NUM,
                out_csv
        };

        cout << "Run test with sorted array." << endl;
        run_testset1(true, test_params1);
        out_csv.close();
    }

    if(test_2) {
        ofstream out_csv_mtls;
        out_csv_mtls.open("parallel_strategies_test_mtls.csv");


        test_params test_params2 = {
                nsizes,
                pb_sizes,
                "speedup parallel_reduce/sequential",
                TESTS_NUM,
                out_csv_mtls
        };

        run_testset2(test_params2);

        out_csv_mtls.close();
    }


    if(test_lcs) {
        ofstream out_csv_lcs;
        out_csv_lcs.open("parallel_strategies_test_lcs.csv");

        test_params test_params3 = {
                nsizes,
                pb_sizes,
                "performance of tiled lcs vs. sequential optim",
                TESTS_NUM,
                out_csv_lcs
        };

        run_testset3(test_params3);

        out_csv_lcs.close();
    }
    delete pb_sizes;
    return 0;
}

int main(int argc, char** argv){
    ofstream out_csv_lcs;
    out_csv_lcs.open("parallel_strategies_test_lcs.csv");


    int nsizes = 2;
    iter_type* pb_sizes = new iter_type[nsizes];
    int start_pow2_size = 10;

    for (int j = 0; j < nsizes; ++j) {
        pb_sizes[j] = (1 << start_pow2_size) + 10000 * j;
    }


    test_params tp = {
            2,
            pb_sizes,
            "performance of tiled lcs vs. sequential optim",
            TESTS_NUM,
            out_csv_lcs
    };

    static task_scheduler_init
            init(task_scheduler_init::deferred);

    init.initialize(LCS_NUM_THREADS, UT_THREAD_DEFAULT_STACK_SIZE);

    tp.out << tp.test_names << endl;
    for (int i = 0; i < tp.nsizes; ++i) {
        long pbsize = tp.pb_sizes[i];
        LCS lcs(pbsize, pbsize);
        lcs.measure_perfs(10);
        LCS_timedata tmd = lcs.get_perfs();
        csvline(tp.out, tmd);
    }

    out_csv_lcs.close();

}