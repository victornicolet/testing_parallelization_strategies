//
// Created by nicolet on 13/09/17.
//

#include "ExamplesTaskBased.h"
#include "Utils.h"
#include "tbb/flow_graph.h"
#include "Stopwatch.h"
#include <tbb/tbb.h>


using namespace tbb::flow;
using namespace tbb;
using namespace std;

#define CHUNK_SIZE 1024

typedef struct {
    iter_type x;
    iter_type y;
    iter_type size;
} square;

typedef struct {
    data_type msqs;
    data_type tsqs;
} msqs_res;

//  Max-top-left square example.
// Three different types of tasks, diagonal, left to diagonal and top to diagonal.

// Sums the element in each line and store them in aux
static void computeLeftElement(data_type** M, data_type* aux, square s) {
    data_type lsum;
    for (iter_type i = s.x; i < s.x + s.size; ++i) {
        lsum = 0;
        for (iter_type j = s.y; j < s.y + s.size; ++j) {
            lsum += M[i][j];
        }
        aux[i-s.x] = lsum;
    }
}

// Sums the element in each column and store them in aux
static void computeTopElement(data_type** M, data_type* aux, square s) {
    data_type lsum;
    for (iter_type j = s.y; j < s.y + s.size; ++j) {
        lsum = 0;
        for (iter_type i = s.x; i < s.x + s.size; ++i) {
            lsum += M[i][j];
        }
        aux[j-s.y] = lsum;
    }
}

// Computes the sums of the borders of top-left squares and store in aux
static void computeDiagElement(data_type ** M, data_type* aux, square s) {
    data_type lsum;

    iter_type b = s.x;
    iter_type e = s.x + s.size;

    for (iter_type i = b; i < e; ++i) {
        lsum = M[i][i];
        for (iter_type j = b; j < i; ++j) {
            lsum += M[j][i] + M[i][j];
        }
        aux[i-b] = lsum;
    }
}

static void joinSum(data_type* auxl, data_type* auxr, iter_type n) {
    for (iter_type i = 0; i < n; ++i) {
        auxl[i] += auxr[i];
    }
}

static msqs_res joinAll(data_type* left, data_type* top, data_type* diag, iter_type n) {
    data_type msqs = 0;
    data_type tsqs = 0;
    for (iter_type i = 0; i < n; ++i) {
        tsqs += left[i] + top[i] + diag[i];
        msqs = std::max(msqs, tsqs);
    }

    return {msqs, tsqs};
}

// The parallel part


struct computeTask {
    data_type **input;
    data_type *aux;
    square s;

    computeTask(data_type** _inp, square _s) {
        input = _inp;
        s = _s;
        aux = new data_type[_s.size];
    }
    void operator()( continue_msg ) const {
        if(s.x == s.y) {
            computeDiagElement(input, aux, s);
        } else if (s.x > s.y) {
            computeLeftElement(input, aux, s);
        } else if (s.x < s.y) {
            computeTopElement(input, aux, s);
        }
    }
};

struct joinTask {
    data_type *auxl;
    data_type *auxr;
    square s;

    void operator()(continue_msg) const {
        joinSum(auxl, auxr, s.size);
    }
};

graph createTaskGraph(data_type **input, iter_type size) {
    int nsquares = (int) (size / CHUNK_SIZE) + 1;
    computeTask** chunk = new computeTask*[nsquares];
    for(int i = 0; i < nsquares; i++) {
        for(int j = 0; j < nsquares; j++) {
            chunk[i] = new computeTask(input,
                                       {i*CHUNK_SIZE,
                                        j*CHUNK_SIZE,
                                        min((i+1) * CHUNK_SIZE, (int) size) - i * CHUNK_SIZE});
        }
    }

}

void taskGraph(data_type** input) {

}

// Parallel reduce style operation
struct MaxTopLeftSquare {
    data_type **A;
    data_type mtops;
    data_type topsum;

    MaxTopLeftSquare(data_type** _in) : mtops(0), topsum(0), A(_in) {}
    MaxTopLeftSquare(MaxTopLeftSquare& mtls, split) : mtops(0), topsum(0), A(mtls.A) {}

    void operator()(const blocked_range<iter_type>& r) {
        for(iter_type i = r.begin(); i != r.end(); i++) {
            topsum += A[i][i];
            for(iter_type j = 0; j < i; j++) {
                topsum += A[i][j] + A[j][i];
            }
            mtops = max(mtops, topsum);
        }
    }

    void join(MaxTopLeftSquare& rhs) {
        topsum += rhs.topsum;
        mtops = max(mtops, topsum + rhs.mtops);
    }
};


double testMaxTopLeftSquareReduction(data_type **in, iter_type n, test_params tp) {
    StopWatch t;
    MaxTopLeftSquare mtls(in);
    double* times = new double[tp.number_per_test];
    for (int i = 0; i < tp.number_per_test; ++i) {
        t.start();
        parallel_reduce(blocked_range<iter_type>(0,n), mtls);
        times[i] = t.stop();
    }

    // Measure sequential
    data_type** A = in;
    data_type mtops = 0;
    data_type topsum = 0;
    t.clear();
    t.start();
    for(iter_type i = 0; i < n; i++) {
        topsum += A[i][i];
        for(iter_type j = 0; j < i; j++) {
            topsum += A[i][j] + A[j][i];
        }
        mtops = max(mtops, topsum);
    }
    double seqtime = t.stop();

    return seqtime / dmean(times, tp.number_per_test);
}