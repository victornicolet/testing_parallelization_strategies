//
// Created by nicolet on 17/01/18.
//

#include <iostream>
using namespace std;

#ifndef PARALLEL_STRATEGIES_TESTING_LARGEST_COMMON_SUBSEQUENCE_H
#define PARALLEL_STRATEGIES_TESTING_LARGEST_COMMON_SUBSEQUENCE_H

#define LCS_NUM_THREADS 16

struct LCS_timedata {
    long pb_size;
    double seq_naive;
    double seq_optim;
    double par_tiled;
    double row_blocks;
    double row_blocks_bis;
};

static void csvline(ostream& out, LCS_timedata tdm){
    out << tdm.pb_size << "," << tdm.seq_naive << "," << tdm.seq_optim << "," << tdm.par_tiled << ","
        << tdm.row_blocks << "," << tdm.row_blocks_bis << endl;
}

class LCS {
    char* X;
    char* Y;
    long X_size;
    long Y_size;

    LCS_timedata perfs;

    long lcs_sequential_naive();
    long lcs_sequential_fast();
    long lcs_parallel_tiled();
    long lcs_parallel_rowblocks();
    long lcs_parallel_rowblocks_constjoin();
    void do_perf_update();
    void print_perfs();

public:
    LCS(long m, long n);
    LCS(char* x, char* y, long m, long n) : X(x),Y(y), X_size(m), Y_size(n) {}
    long longest_common_subsequence(int strategy);
    void measure_perfs(int n);
    LCS_timedata get_perfs() { return perfs; }
};


#endif //PARALLEL_STRATEGIES_TESTING_LARGEST_COMMON_SUBSEQUENCE_H
