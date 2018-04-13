//
// Created by victorn on 13/04/18.
//

#ifndef PARALLEL_STRATEGIES_TESTING_MAXTOPRIGHTRECTANGLE_H
#define PARALLEL_STRATEGIES_TESTING_MAXTOPRIGHTRECTANGLE_H

using namespace std;

#include "Utils.h"


#define mTRR_NUM_THREADS 32
#define DATA_RANGE 10

typedef struct etimes {
    int threads;
    double seq;
    double split;
    double fused;
} etimes;

class MaxTopRightRectangle {
    int** A;
    int* col;
    long n;
    long m;
    std::list<etimes> experiments;

public:
    int mtrr;

    MaxTopRightRectangle(long m, long n): mtrr(0), m(m), n(n){
        cout << "Initialization m = " << m << " n = " << n << endl;
        A = new int*[n];
        for(long i = 0; i < n; i++){
            A[i] = new int[m];
            cout << "\r" << (i * 100 / n) << "%";
            for(long j = 0; j < m; j++){
                A[i][j] = (rand() % DATA_RANGE) - (DATA_RANGE / 2);
            }
        }
        cout << endl;
        col = new int[m];

        for(long j = 0; j < m; j++){
            col[0] = 0;
        }
    }

    int run_exp(int nthreads, int seqruns, int splitruns, int fusedruns);
    int mtrr_seq();
    int mtrr_parallel(bool fused);
    int output(ofstream& out);
};


#endif //PARALLEL_STRATEGIES_TESTING_MAXTOPRIGHTRECTANGLE_H
