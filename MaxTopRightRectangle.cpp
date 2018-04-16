//
// Created by victorn on 13/04/18.
//

#include "MaxTopRightRectangle.h"

#include "Stopwatch.h"

using namespace std;
/* All the kernels for the max-top-right rectangle. We need an auxiliary for the
 * parallel case so we have different kernels for the sequential case and the
 * parallel case. For the parallel case, we also have two different kernels. One
 * is forced to use only one pass, the other one would be the original sequential
 * kernel + another loop in the other direction to do the max-terminal-sums.
*/

void _mtrr_seq_kernel_ (int** A, int* col, int& mtrr, long m, long i0, long iend){
    int rmr = 0;

    long i,j;
    for(i = i0; i < iend; i++){
        rmr = 0;
        for(j = 0; j < m; j++){
            col[j] += A[i][j];
            rmr = max(col[j] + rmr, 0);
        }
        mtrr = max(rmr, mtrr);
    }
}

void _mtrr_par_split_kernel(int** A, int* col, int* mtsj, int& mtrr, long m, long i0, long iend){
    int rmr = 0;
    int col_sums;

    int** _a = A;

    for(long i = i0; i < iend; i++){
        rmr = 0;

//         Rightwards original loop
        for(long j = 0; j < m; j++){
            col[j] += _a[i][j];
            rmr = max(col[j] + rmr, 0);
        }

#ifndef MTRR_LEFTWARDS_AUX
//        Leftwards auxiliary loop
        col_sums = 0;
        for(long j = m-1; j >= 0; j--){
            col_sums += col[j];
            mtsj[j] = max(col_sums, mtsj[j]);
        }
#else
// Experiment with both loops in the same direction.
        col_sums = 0;
        for(long j = 0; j < m; j++){
            col_sums += col[j];
            mtsj[j] = max(col_sums, mtsj[j]);
        }
#endif

        mtrr = max(rmr, mtrr);
    }
}

void _mtrr_par_fused_kernel(int** A, int* col, int* mtsj, int& mtrr, long m, long i0, long iend){
    int rmr = 0;
    int sumcol = 0;

    long i;
    for(i = i0; i < iend; i++){
        rmr = 0;
//        Leftwards auxiliary loop + fused original loop
        sumcol = 0;
        for(long j = m-1; j >= 0; j--){
            col[j] += A[i][j];
            sumcol += col[j];
            rmr = max(sumcol, rmr);
            mtsj[j] = max(sumcol, mtsj[j]);
        }

        mtrr = max(rmr, mtrr);
    }
}


void _mtrr_join_ (int* col_l, int* col_r, int* mtsj_l, int* mtsj_r, int& mtrr_l,
                  int& mtrr_r, long m) {
    long j = 0;
    int tersum = 0;
    mtrr_l = 0;
    for(j = m-1; j>=0; j--){
        tersum += col_l[j];
        mtsj_l[j] = max(mtsj_r[j] + tersum, mtsj_l[j]);
        col_l[j] += col_r[j];
        mtrr_l = max(mtrr_l, mtsj_l[j]);
    }
}


struct MtrrReduction {
    int** A;
    int* col;
    int* mtsj;
    int mtrr;
    long m;
    bool fused;

    MtrrReduction(int** A, int* col, int* mtsj, long m, int mtrr, bool fused) :
            A(A), col(col), mtsj(mtsj), m(m), mtrr(mtrr), fused(fused) {}

    MtrrReduction(MtrrReduction& s, split) {
        A = s.A;
        col = s.col;
        mtsj = s.mtsj;
        mtrr = s.mtrr;
        m = s.m;
        fused = s.fused;
    }

    void operator()(const blocked_range<long>& r){
        if(fused)
            _mtrr_par_fused_kernel(A, col, mtsj, mtrr, m, r.begin(), r.end());
        else
            _mtrr_par_split_kernel(A, col, mtsj, mtrr, m, r.begin(), r.end());
    }

    void join(MtrrReduction& rhs) {
        _mtrr_join_(col, rhs.col, mtsj, rhs.mtsj, mtrr, rhs.mtrr, m);
    }
};


int MaxTopRightRectangle::mtrr_parallel(bool fused) {
    auto mtsj = new int[m];
    for(long j = 0; j < m; j++){
        mtsj[j] = 0;
        col[j] = 0;
    }

    MtrrReduction mtrr_tbb(A, col, mtsj, m, mtrr, fused);
    parallel_reduce(blocked_range<long>(0,n), mtrr_tbb);
    mtrr = mtrr_tbb.mtrr;

    return mtrr;
}

int MaxTopRightRectangle::mtrr_seq() {
    for(long j = 0; j < m; j++){
        col[j] = 0;
    }
    _mtrr_seq_kernel_(A, col, mtrr, m, 0, n);
}

int MaxTopRightRectangle::run_exp(int nthreads, int seqruns, int splitruns, int fusedruns) {
    int i = 0;
    StopWatch t;
//    Sequential
    cout << "Sequential experiment ...." << endl;
    auto times = new double[seqruns];
    for(i = 0; i < seqruns; i++){
        t.start();
        mtrr = this->mtrr_seq();
        times[i] = t.stop();
    }
//    Parallel, but loops are not fused
    cout << "Parallel split experiment ...." << endl;
    auto ptimes = new double[splitruns];
    for(i = 0; i < splitruns; i++){
        t.start();
        mtrr = this->mtrr_parallel(false);
        ptimes[i] = t.stop();
    }

// Parallel with fused loops.
    cout << "Parallel fused experiment ...." << endl;
    auto p0times = new double[fusedruns];
    for(i = 0; i < fusedruns; i++){
        t.start();
        mtrr = this->mtrr_parallel(true);
        p0times[i] = t.stop();
    }

    etimes experiment_results{ nthreads,
                           dmean(times, seqruns),
                           dmean(ptimes, splitruns),
                           dmean(p0times, fusedruns)};

    experiments.push_back(experiment_results);
}

int MaxTopRightRectangle::output(ofstream& out) {
    for(auto const& i: experiments){
        out << n << "," << m << "," << i.threads << "," << i.seq << "," << i.split << "," << i.fused << endl;
    }
}