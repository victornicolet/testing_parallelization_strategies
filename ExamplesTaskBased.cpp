//
// Created by nicolet on 13/09/17.
//

#include "ExamplesTaskBased.h"
#include "Utils.h"
#include "tbb/flow_graph.h"
#include "Stopwatch.h"
#include <tbb/tbb.h>

#define CHUNK_SIZE 1024

using namespace tbb::flow;
using namespace tbb;
using namespace std;

typedef continue_node< continue_msg > node_t;
typedef const continue_msg & msg_t;

typedef struct flow_graph {
    broadcast_node<continue_msg> start_node;
    graph* g;
} flow_graph;

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

static msqs_res join_msqs(msqs_res tl, msqs_res cur) {
    return {tl.tsqs + cur.tsqs, max(tl.msqs, tl.tsqs + cur.msqs)};
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

graph createTaskGraphFullPar(data_type **input, iter_type size) {
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

// Other implementation: graph task, but not fully parallel

struct fullComputeTask {
    data_type **input;
    data_type *aux;
    square s;
    msqs_res m;
    fullComputeTask* leftParent;
    fullComputeTask* topParent;
    fullComputeTask* topLeftParent;

    fullComputeTask(data_type** _inp, square _s) {
        m = {0,0};
        input = _inp;
        s = _s;
        aux = new data_type[_s.size];
    }

    fullComputeTask(fullComputeTask* p, square _s) {
        m = {0,0};
        input = p->input;
        s = _s;
        aux = new data_type[_s.size];

        if(s.x == s.y) {
            topLeftParent = p;
        } else if (s.x > s.y) {
            leftParent = p;
        } else if (s.x < s.y) {
            topParent = p;
        }
    }

    fullComputeTask(fullComputeTask* pl, fullComputeTask* ptl, fullComputeTask* pt, square _s) {
        m = {0,0};
        input = pl->input;
        s = _s;
        aux = new data_type[_s.size];
        topLeftParent = ptl;
        leftParent = pl;
        topParent = pt;
    }

    void operator()( continue_msg ) {
        if(s.x == s.y) {
            computeDiagElement(input, aux, s);
            m = join_msqs(topLeftParent->m,
                                   joinAll(leftParent->aux, topParent->aux, aux, s.size));
        } else if (s.x > s.y) {
            computeLeftElement(input, aux, s);
            joinSum(leftParent->aux, aux, s.size);
        } else if (s.x < s.y) {
            computeTopElement(input, aux, s);
            joinSum(topParent->aux, aux, s.size);
        }
    }
};

flow_graph createTaskGraphPipelined(data_type **input, iter_type size) {
    int nsquares = (int) (size / CHUNK_SIZE) + 1;
    graph g;



    fullComputeTask*** chunk = new fullComputeTask**[nsquares];
    node_t*** nodes = new node_t**[nsquares];
    // Create one empty parent task
    fullComputeTask* t0 = new fullComputeTask(input, {0,0,0});
    broadcast_node< continue_msg > parent(g);

    if (nsquares < 10) {
        cout << "Problem is too small! Won't do anything." << endl;
        return {parent, &g};
    }

    // Create the first tasks (top line, first diagonal, left column)
    for(int i = 0; i < nsquares; i++) {
        chunk[i] = new fullComputeTask*[nsquares];
        nodes[i] = new node_t*[nsquares];
    }

    // The first task
    chunk[0][0] = new fullComputeTask(t0,t0,t0,{0,0,CHUNK_SIZE});
    nodes[0][0] = new node_t(g, chunk[0][0]);
    make_edge(parent, *nodes[0][0]);

    // Top line and column
    for(int i = 1; i < nsquares; i++) {
        chunk[0][i] = new fullComputeTask(t0, {i*CHUNK_SIZE, 0,
                                                min((i+1) * CHUNK_SIZE, (int) size) - i * CHUNK_SIZE});
        nodes[0][i] = new node_t(g, chunk[0][i]);
        make_edge(parent, *nodes[0][i]);

        chunk[i][0] = new fullComputeTask(t0, {0, i*CHUNK_SIZE,
                                               min((i+1) * CHUNK_SIZE, (int) size) - i * CHUNK_SIZE});
        make_edge(parent, *nodes[i][0]);
    }

    // All dependent tasks
    for(int i = 1; i < nsquares; i++) {
        chunk[i] = new fullComputeTask*[nsquares];
        for(int j = 1; j < nsquares; j++) {
            if(i == j) {
                // This is a diagonal task
                chunk[i][j] = new fullComputeTask(chunk[i][j-1],
                                                  chunk[i-1][j-1],
                                                  chunk[i-1][j],
                                                  {i * CHUNK_SIZE,
                                                   j * CHUNK_SIZE,
                                                   min((i + 1) * CHUNK_SIZE, (int) size) - i * CHUNK_SIZE});

            } else if (i < j) {
                // This is a line task
                chunk[i][j] = new fullComputeTask(chunk[i][j-1],
                                                  {i * CHUNK_SIZE,
                                                   j * CHUNK_SIZE,
                                                   min((i + 1) * CHUNK_SIZE, (int) size) - i * CHUNK_SIZE});
            } else {
                // This is a column task
                chunk[i][j] = new fullComputeTask(chunk[i-1][j],
                                                  {i * CHUNK_SIZE,
                                                   j * CHUNK_SIZE,
                                                   min((i + 1) * CHUNK_SIZE, (int) size) - i * CHUNK_SIZE});
            }

            // Create the nde corresponding to the task
            nodes[i][j] = new node_t(g, chunk[i][j]);

            // Create the dependency edges on the graph
            if(i==j) {
                make_edge(*nodes[i-1][j], *nodes[i][j]);
                make_edge(*nodes[i][j-1], *nodes[i][j]);
                make_edge(*nodes[i-1][j-1], *nodes[i][j]);
            } else if (i < j) {
                make_edge(*nodes[i][j-1], *nodes[i][j]);
            } else {
                make_edge(*nodes[i][j-1], *nodes[i][j]);
            }

        }
    }

    return {parent, &g};
}

result_data testMaxTopLeftSquareTaskPipelined (data_type** in, iter_type n, test_params tp){
    StopWatch t;

    double* times = new double[tp.number_per_test];
    for (int i = 0; i < tp.number_per_test; ++i) {
        flow_graph fg = createTaskGraphPipelined(in, n);
        t.start();
        fg.start_node.try_put(continue_msg());
        fg.g->wait_for_all();
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
    double st1par_time = dmean(times, tp.number_per_test);
    cout << "Speedup: " << seqtime / st1par_time;

    return {seqtime, st1par_time, 0.0, n, "maxTopLeftSquare with pipelined tasks"};
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


result_data testMaxTopLeftSquareReduction(data_type **in, iter_type n, test_params tp) {
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
    double st1par_time = dmean(times, tp.number_per_test);
    cout << "Speedup: " << seqtime / st1par_time;

    return {seqtime, st1par_time, 0.0, n, "maxTopLeftSquare"};
}