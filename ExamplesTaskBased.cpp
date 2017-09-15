//
// Created by nicolet on 13/09/17.
//

#include "ExamplesTaskBased.h"
#include "Stopwatch.h"

#define CHUNK_SIZE 256

using namespace tbb::flow;
using namespace tbb;
using namespace std;

typedef continue_node< continue_msg > node_t;

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
struct taskRes {
public:
    data_type sum;
    data_type maxs;

    taskRes(data_type sum, data_type m) : sum(sum), maxs(m) {}
};


static void join_msqs(taskRes* tl, taskRes* cur, taskRes* m) {
    m->sum = tl->sum + cur->sum;
    m->maxs =  max(tl->sum, tl->sum + cur->maxs);
}

struct fullComputeTask {
    data_type **input;
    data_type *aux;
    square s;
    taskRes *m;
    fullComputeTask* leftParent;
    fullComputeTask* topParent;
    fullComputeTask* topLeftParent;

    fullComputeTask(data_type** _inp, square _s) {
        m = new taskRes(0,0);
        input = _inp;
        s = _s;
        aux = new data_type[_s.size];
    }

    fullComputeTask(fullComputeTask* p, square _s) {
        m = new taskRes(0,0);
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
        m = new taskRes(0,0);
        input = pl->input;
        s = _s;
        aux = new data_type[_s.size];
        topLeftParent = ptl;
        leftParent = pl;
        topParent = pt;
    }

    void operator()( continue_msg ) const {
        if(s.x == s.y) {
            computeDiagElement(input, aux, s);
            msqs_res _msqs = joinAll(leftParent->aux, topParent->aux, aux, s.size);
            join_msqs(topLeftParent->m, new taskRes(_msqs.tsqs, _msqs.msqs),m);
        } else if (s.x > s.y) {
            computeLeftElement(input, aux, s);
            joinSum(leftParent->aux, aux, s.size);
        } else if (s.x < s.y) {
            computeTopElement(input, aux, s);
            joinSum(topParent->aux, aux, s.size);
        }
    }
};

void populateTaskGraphPipelined(graph* g, node_t parent,
                                data_type **input, iter_type size) {
    int nsquares = (int) (size / CHUNK_SIZE) + 1;

    fullComputeTask*** chunk = new fullComputeTask**[nsquares];
    node_t*** nodes = new node_t**[nsquares];
    // Create one empty parent task
    fullComputeTask* t0 = new fullComputeTask(input, {0,0,0});

    if (nsquares < 10) {
        cout << "Problem is too small! Won't do anything." << endl;

    }

    // Create the first tasks (top line, first diagonal, left column)
    for(int i = 0; i < nsquares; i++) {
        chunk[i] = new fullComputeTask*[nsquares];
        nodes[i] = new node_t*[nsquares];
    }

    // The first task
    chunk[0][0] = new fullComputeTask(t0,t0,t0,{0,0,CHUNK_SIZE});
    nodes[0][0] = new node_t(*g, *chunk[0][0]);
    make_edge(parent, *nodes[0][0]);

    // Top line and column
    for(int i = 1; i < nsquares; i++) {
        chunk[0][i] = new fullComputeTask(t0, {i*CHUNK_SIZE, 0,
                                                min((i+1) * CHUNK_SIZE, (int) size) - i * CHUNK_SIZE});
        nodes[0][i] = new node_t(*g, *chunk[0][i]);
        make_edge(parent, *nodes[0][i]);

        chunk[i][0] = new fullComputeTask(t0, {0, i*CHUNK_SIZE,
                                               min((i+1) * CHUNK_SIZE, (int) size) - i * CHUNK_SIZE});
        nodes[0][i] = new node_t(*g, *chunk[i][0]);
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
            nodes[i][j] = new node_t(*g, *chunk[i][j]);

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
}

result_data testMaxTopLeftSquareTaskPipelined (data_type** in, iter_type n, test_params tp){
    StopWatch t;

    double* times = new double[tp.number_per_test];
    graph g;
    node_t start_node(g, [](continue_msg){} );
    populateTaskGraphPipelined(&g, start_node, in, n);


    for (int i = 0; i < tp.number_per_test; ++i) {
        t.start();
        start_node.try_put(continue_msg());
        g.wait_for_all();
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
    cout << "Speedup: " << seqtime / st1par_time << endl;

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
    cout << "Speedup reduction naive/ seq naive: " << seqtime / st1par_time << endl;

    return {seqtime, st1par_time, 0.0, n, "maxTopLeftSquare"};
}



// Other parallel solution using scans and only one outer loop parallelized
// Possibly could use one row for the colsums and avoid copying but for now
// don't care and jsut copy along.
struct MtlsMultiscan {
    data_type **A;
    iter_type b,e, n;
    // Linear size temp storage
    data_type *colsums;
    data_type* rowsums;

    MtlsMultiscan(data_type** _in, iter_type n ) : n(n), b(-1), e(-1), A(_in) {
        colsums = new data_type[n];
        rowsums = new data_type[n];
    }

    MtlsMultiscan(MtlsMultiscan& mtls, split) : A(mtls.A) {
        b= -1; e= -1; n = mtls.n;
        colsums = new data_type[n];
        rowsums = new data_type[n];
    }

    void operator()(const blocked_range<iter_type>& r) {
        data_type rsum = 0;
        for(iter_type i = r.begin(); i != r.end(); i++) {
            rsum = 0;
            for(iter_type j = 0; j < i; j++) {
                rsum += A[i][j];
            }
            for(iter_type j = i + 1; j < n; j++) {
                colsums[j] += A[i][j];
            }
            rowsums[i] = rsum;
        }
        b = r.begin();
        e = r.end();
    }

    void join(MtlsMultiscan& rhs) {

        // Cells of the righthand's side just need to be copied into the leftside's array
        for (iter_type i = rhs.b; i < rhs.e; i++) {
            rowsums[i] = rhs.rowsums[i];
        }
        // Add the sums accumulated so far in the array.
        for (iter_type j = rhs.b; j < n; j++) {
            colsums[j] += rhs.colsums[j];
        }
        // Update the borders
        e = rhs.e;
    }
};

struct ReduceMultiScanProd {
    iter_type n;
    data_type *rowsums;
    data_type *colsums;
    data_type sum;
    data_type mtls;

    ReduceMultiScanProd(data_type* rs, data_type* cs, iter_type n) :
            rowsums(rs), colsums(cs), n(n), sum(0), mtls(0) {}

    ReduceMultiScanProd(ReduceMultiScanProd& rmp, split) : sum(0), mtls(0) {
        rowsums = rmp.rowsums;
        colsums = rmp.colsums;
    }

    void operator()(blocked_range<iter_type>& r) {
        for (iter_type i = r.begin(); i < r.end(); ++i) {
            sum += rowsums[i] + colsums[i];
            mtls = max(mtls, sum);
        }
    }

    void join(ReduceMultiScanProd& rhs) {
        mtls = max(mtls, sum + rhs.mtls);
        sum += rhs.sum;
    }

};

data_type mtlsMultiscanSequential(data_type** M, iter_type n){
    data_type mtls, sum, aux_sum;
    data_type* aux_rows, *aux_cols;

    // Fill in the aux_rows and aux_cols with the partial sums
    aux_rows = new data_type[n];
    aux_cols = new data_type[n];
    for(iter_type i = 0 ; i < n; i ++) {
        aux_sum = M[i][i];
        for(iter_type j = 0; j < i; j++) {
            aux_sum += M[i][j];
        }
        for(iter_type j = i + 1; j < n; j++) {
            aux_cols[j] += M[i][j];
        }
        aux_rows[i] = aux_sum;
    }

    // Compute the mtls using the sums stored in the aux arrays.
    sum = 0;
    mtls = 0;
    for(iter_type i = 0; i <n; i++) {
        sum += aux_rows[i] + aux_cols[i];
        mtls = max(mtls, sum);
    }

    return mtls;
}

result_data testMtlsMultiscan(data_type** M, iter_type n, test_params tp) {
    StopWatch t;
    MtlsMultiscan mtls(M, n);
    double* times = new double[tp.number_per_test];
    double mtls_par_opt = 0;
    for (int i = 0; i < tp.number_per_test; ++i) {
        t.start();
        // First perform the column and row sums
        parallel_reduce(blocked_range<iter_type>(0,n), mtls);
        // Reduce them like for the mps
        ReduceMultiScanProd rmsp(mtls.rowsums, mtls.colsums,n);
        parallel_reduce(blocked_range<iter_type>(0,n), rmsp);
        times[i] = t.stop();
        mtls_par_opt = rmsp.mtls;
    }

    // Measure sequential
    data_type** A = M;
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


    // Measure "optimized" sequential.
    t.clear();
    t.start();
    double mtls_seq_opt = mtlsMultiscanSequential(M, n);
    double mtls_seq_opt_time = t.stop();

    cout << "Speedup multiscan parallel / seq_naive: " << seqtime / st1par_time << endl;
    cout << "Speedup seq_opt/seq_naive : " << seqtime/mtls_seq_opt_time << endl;
    cout << "Test mtls_seq_opt == mtls_par_opt ? " << (mtls_par_opt == mtls_seq_opt) << endl;

    return {seqtime, st1par_time, mtls_seq_opt_time, n, "MtlsMultiscan"};
}