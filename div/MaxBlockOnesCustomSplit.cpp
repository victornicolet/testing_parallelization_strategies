//
// Created by victorn on 11/09/18.
//

#include <iostream>
#include <random>
#include <climits>
#include "omp.h"

// Parallelism parameters
#define CHUNKSIZE (1024 * 8 * 4)
#define UB 255.0
#define LB 0.0
#define COND 127.5

using namespace std;

struct LongChunk {
    long start;
    long end;
    long end_prefix;
};


long mymax(long r, long n) {
    // r is the already reduced value
    // n is the new value
    long m;
    if (n > r) {
        m = n;
    } else {
        m = r;
    }
    return m;
}

double f(double x) {
    return sqrt(x);
}

long _longest_block(double *A, long n, int nthreads) {

    omp_set_dynamic(0);
    omp_set_num_threads(nthreads);


    long num_chunks = n / CHUNKSIZE + 1;


    LongChunk chunks[num_chunks];
    long results[num_chunks];

    double t0 = omp_get_wtime();

#pragma omp parallel for
    for (long cno = 0; cno < num_chunks; cno++) {
        chunks[cno].start = min(cno * CHUNKSIZE, n);
        long end = min((cno + 1) * CHUNKSIZE, n);
        long end0 = end;
        bool prefix = A[end0 + 1] < COND;
        while(prefix && end0 < max(n, (cno + 2)* CHUNKSIZE)) {
            prefix = prefix && A[end0] < COND;
            end0++;
        }

        chunks[cno].end = end;
        chunks[cno].end_prefix = end0;
    }

#pragma omp parallel for
    for (long cno = 0; cno < num_chunks; cno++) {
        long st = chunks[cno].start;
        long en = chunks[cno].end_prefix;
        long ml = 0;
        long cl = 0;
        for (long i = st; i < en; i++) {
            cl = (f(A[i]) < COND) ? cl + 1 : cl;
            ml = max(ml, cl);
        }

        results[cno] = ml;
    }


#pragma omp declare reduction \
  (rwz:long:omp_out=mymax(omp_out, omp_in)) \
  initializer(omp_priv=0)

    long m = 0;


#pragma omp parallel for reduction(rwz:m)
    for (long cno = 0; cno < num_chunks; cno++)
        m = mymax(m, results[cno]);

    double t_final = omp_get_wtime();

    cout << "Elapsed:" << t_final - t0 << endl;


    long sum_pl = 0;
    for (long cno = 0; cno < num_chunks; cno++)
        sum_pl += chunks[cno].end_prefix - chunks[cno].end;

    float avg_prefix = (float) sum_pl / num_chunks;

    float dev = 0.0;

    for (long cno = 0; cno < num_chunks; cno++)
        dev += ((float)(chunks[cno].end_prefix - chunks[cno].end)) - avg_prefix;

    cout << "Average prefix length: " << avg_prefix << endl;
    cout << "Deviance: " << dev << endl;

    return m;
}


int main(int argc, char** argv) {
    long n = 1 << atoi(argv[1]);
    int t = atoi(argv[2]);

    double* A = (double*)malloc(sizeof(double) * n);

//     Array creation and initialization
    double lb = LB;
    double ub = UB;
    uniform_real_distribution<double> uniformRealDistribution(lb, ub);

    default_random_engine re;

#pragma omp parallel for
    for (long j = 0; j < n; j++) {
        A[j] = uniformRealDistribution(re);
    }



    long ml = _longest_block(A, n, t);

    return 0;
}