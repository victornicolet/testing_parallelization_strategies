//
// Created by nicolet on 17/01/18.
//
#include <math.h>
#include <cstdlib>
#include "tbb/tbb.h"
#include <cstdio>
#include "longest_common_subsequence.h"
#include <time.h>
#include "Stopwatch.h"


#define TILE_SIZE 1024

using namespace tbb;


static char GenerateRandomChar() {
    return (char) 'A' + rand()%24;
}

static long max(long a, long b){
    return a > b ? a : b;
}

static long min(long a, long b){
    return a > b ? b : a;
}

LCS::LCS(long m, long n) {
    X_size = m;
    Y_size = n;

    perfs.pb_size = max(m,n);
    perfs.seq_naive = 0.0;
    perfs.seq_optim = 0.0;
    perfs.par_tiled = 0.0;

    srand(time(NULL));
//    Initialize the data randomly
    X = new char[m];
    Y = new char[n];
    for(long i = 0; i< m; i++){
        X[i] = GenerateRandomChar();
    }
    for(long j = 0; j < n; j++){
        Y[j] = GenerateRandomChar();
    }
}

static long** alloc_aux_matrix(long m, long n){
    long** C = new long*[m+1];
    long i,j;

    for(i = 0; i <= m; i++){
        C[i] = new long[n+1];
        C[i][0] = 0;
    }
    for(i = 1; i <= n; i++){
        C[0][i] = 0;
    }

    return C;
}

static void delete_aux_matrix(long** C, long m) {
    long i;
    for(i = 0; i <= m; i ++) {
        delete [] C[i];
    }
    delete [] C;
    return;
}

inline void _lcs_kernel_reg_ (long bx, long ex, long by, long ey, long** C, char *X, char* Y){
    long i,j;
    for(i = bx; i < ex; i++){
        for(j = by; j < ey; j++) {
            if (X[i - 1] == Y[j - 1]) {
                C[i][j] = C[i - 1][j - 1] + 1;
            } else {
                C[i][j] = max(C[i - 1][j], C[i][j - 1]);
            }
        }
    }
}
long LCS::lcs_sequential_naive(){
    long lcs ;
    long m = X_size;
    long n = Y_size;

    long ** C = alloc_aux_matrix(m,n);
// Just execute the kernel on the whole matrix
    _lcs_kernel_reg_(1, m+1, 1, n+1, C, X, Y);

    lcs = C[m][n];
    delete_aux_matrix(C,m);
    return lcs;
}


// An amelioration of the sequential version using an array as auxiliary.
long LCS::lcs_sequential_fast(){
    long lcs = 0L;
    long i,j;
    long r0 = 0;
    long r1 = 0;
    long m = X_size;
    long n = Y_size;
//    Allocate the temporary matrix
    long* C = new long[n+1];
    for(i = 0; i <= n; i++)
        C[i] = 0;

    for(i = 1; i <= m; i++) {
        r1 = 0;
        for(j = 1; j <= n; j++) {
            r0 = C[j];
            if(X[i-1] == Y[j-1])
                C[j] = r1 + 1 ;
            else
                C[j] = max(C[j],C[j-1]);
            r1 = r0;
        }
    }
    lcs = C[n];
    return lcs;
}


// A regular tile kernel


// A parallel tiled version of the sequential algorithm using a matrix for dynamic programming
long LCS::lcs_parallel_tiled(){
    long lcs = 1;
    long **C;
    long i, j;
    long diagnum;
    long cellnum;
    long Ndiags;
    long m = X_size;
    long n = Y_size;

    Ndiags = (n / TILE_SIZE) + (m / TILE_SIZE) + 1;
    long xdiags = (m / TILE_SIZE) + 1;
    long ydiags = (n / TILE_SIZE) + 1;
#ifdef DEBUG_TILES
    printf("m = %ld,  n = %ld\nNdiags = %ld\n", m,n, Ndiags);
#endif
    // ALlocate the auxiliary matrix C
    //    Allocate the temporary matrix
    C = alloc_aux_matrix(m,n);
//    Execute the tiles
    for(diagnum = 0; diagnum <= (Ndiags - 1) * 2; diagnum++) {

        long first_cell = max(diagnum - xdiags - 1, 0);
        long last_cell = min(ydiags + 1, diagnum);
#ifdef DEBUG_TILES
        printf("Diagonal %ld, first cell %ld, last cell %ld\n", diagnum, first_cell, last_cell);
#endif

#pragma parallel for num_threads(LCS_NUM_THREADS)
        for(cellnum = first_cell; cellnum <= last_cell; cellnum++ ) {
//            Tile kernel
            long bx = max(1,(diagnum - cellnum)*TILE_SIZE );
            long ex = min(m + 1, (diagnum - cellnum + 1) * TILE_SIZE);
            long by = max(1, cellnum*TILE_SIZE);
            long ey = min(n + 1, (cellnum + 1) * TILE_SIZE);
#ifdef DEBUG_TILES
            printf("Tile (%ld, %ld) - (%ld, %ld)\n", bx, by, ex, ey);
#endif
            _lcs_kernel_reg_(bx, ex, by, ey, C,X,Y);
        }
    }
    lcs = C[m][n];
    delete_aux_matrix(C,m);
    return lcs;
}


// The optimal? version

void _lcs_super_kernel_linear_ (long start, long end, long n,
                                long* C, long* Cs, long* Cp, bool* cnd,
                                char* X, char* Y){
    long i,j;
    long r0, r1, rs0, rs1;
    bool b;
    C[0] = 0;
    for (i = start; i <= end; ++i) {
        r1 = 0;
        rs1 = 0;
        for(j = 1; j < n + 1; j++) {
            b = X[i-1] == Y[j-1];
//            LCS
            r0 = C[j];
            if(b)
                C[j] = r1 + 1 ;
            else
                C[j] = max(C[j],C[j-1]);
            r1 = r0;
//            Suffix
            rs0 = Cs[j];

            if(b)
                Cs[j] = rs1 + 1;
            else
                Cs[j] = 0;

            rs1 = rs0;
//            Prefix
            if(j - i >= 0){
                cnd[j-i] = cnd[j-i] && b;
                if(cnd[j-i])
                    Cp[j-i] += 1;
                else
                    Cp[j-i] += 0;
            }
        }
    }
    return;
}

struct _lcs {
    long n;
    long* Cs;
    long* Cp;
    long* C;
    long lcs;
};

_lcs _lcs_join_(_lcs left, _lcs right){
    long j;
    long Csp_max = 0;
    for(j = 0; j <= right.n; j++){
        Csp_max = max(Csp_max, left.Cs[j] + right.Cp[j]);
    }
    return {left.n, right.Cs, left.Cp, left.C, max(max(left.lcs, right.lcs), Csp_max)};
}

struct LcsReduction {
    char *X;
    char *Y;
    long n, m;
    bool* cnd;
    _lcs lcs;

    LcsReduction(char* _X, char* _Y, long m, long n) : X(_X), Y(_Y), m(m), n(n) {
        long *C = new long[n+1];
        long *Cp = new long[n+1];
        long *Cs = new long[n+1];
        cnd = new bool[n+1];

        for(long i = 0; i < n +1; i++) {
            C[i] = 0;
            Cp[i] = 0;
            Cs[i] = 0;
            cnd[i] = true;
        }
        lcs = {n, Cs, Cp, C, 0};
    }

    LcsReduction(LcsReduction& s, split): X(s.X), Y(s.Y), m(s.m), n(s.n) {
        long _n = s.n;
        long *C = new long[_n+1];
        long *Cp = new long[_n+1];
        long *Cs = new long[_n+1];
        cnd = new bool[n+1];

        for(long i = 0; i < _n +1; i++) {
            C[i] = 0;
            Cp[i] = 0;
            Cs[i] = 0;
            cnd[i] = true;
        }
        lcs = {_n, Cs, Cp, C, 0};

    }

    void operator()(const blocked_range<long>& r){
        _lcs_super_kernel_linear_(r.begin(), r.end(), n, lcs.C, lcs.Cs, lcs.Cp, cnd, X, Y);
        lcs.lcs = lcs.C[n];
    }

    void join(LcsReduction& rhs) {
        lcs = _lcs_join_(lcs, rhs.lcs);
        for(long i = 0; i < n + 1; i++){
            cnd[i] = true;
            lcs.C[i] = 0;
        }
    }
};



struct LcsReductionConstJoin {
    char *X;
    char *Y;
    long n, m;
    bool* cnd;
    _lcs lcs;

    LcsReductionConstJoin(char* _X, char* _Y, long m, long n) : X(_X), Y(_Y), m(m), n(n) {
        long *C = new long[n+1];
        long *Cp = new long[n+1];
        long *Cs = new long[n+1];
        cnd = new bool[n+1];

        for(long i = 0; i < n +1; i++) {
            C[i] = 0;
            Cp[i] = 0;
            Cs[i] = 0;
            cnd[i] = true;
        }
        lcs = {n, Cs, Cp, C, 0};
    }

    LcsReductionConstJoin(LcsReductionConstJoin& s, split): X(s.X), Y(s.Y), m(s.m), n(s.n) {
        long _n = s.n;
        long *C = new long[_n+1];
        long *Cp = new long[_n+1];
        long *Cs = new long[_n+1];
        cnd = new bool[n+1];

        for(long i = 0; i < _n +1; i++) {
            C[i] = 0;
            Cp[i] = 0;
            Cs[i] = 0;
            cnd[i] = true;
        }
        lcs = {_n, Cs, Cp, C, 0};

    }

    void operator()(const blocked_range<long>& r){
        _lcs_super_kernel_linear_(r.begin(), r.end(), n, lcs.C, lcs.Cs, lcs.Cp, cnd, X, Y);
        lcs.lcs = lcs.C[n];
    }

    void join(LcsReductionConstJoin& rhs) {
        lcs = _lcs_join_(lcs, rhs.lcs);
        for(long i = 0; i < n + 1; i++){
            cnd[i] = true;
            lcs.C[i] = 0;



        }
    }
};


long LCS::lcs_parallel_rowblocks() {
    LcsReduction lcsr(X,Y,X_size,Y_size);
    parallel_reduce(blocked_range<long>(0,X_size), lcsr);
    return lcsr.lcs.lcs;
}

long LCS::lcs_parallel_rowblocks_constjoin() {
    LcsReductionConstJoin lcsr(X,Y,X_size,Y_size);
    parallel_reduce(blocked_range<long>(0,X_size), lcsr);
    return lcsr.lcs.lcs;
}

long LCS::longest_common_subsequence(int strategy){
//    Strategy 1 : sequential naïve implementation
    switch (strategy) {
        case 1:
            return lcs_sequential_naive();
        case 2:
            return lcs_sequential_fast();
        case 3:
            return lcs_parallel_tiled();
        case 4:
            return lcs_parallel_rowblocks();
        case 5:
            return lcs_parallel_rowblocks_constjoin();
        default:
            return 0L;
    }
}

void LCS::do_perf_update() {
    StopWatch t;
    double t1,t2,t3, t4, t4bis;

    t.start();
    lcs_sequential_naive();
    t1 = t.stop();

    t.start();
    lcs_sequential_fast();
    t2 = t.stop();

    t.start();
    lcs_parallel_tiled();
    t3 = t.stop();

    t.start();
    lcs_parallel_rowblocks();
    t4 = t.stop();

    t.start();
    lcs_parallel_rowblocks();
    t4bis = t.stop();

    perfs.seq_naive = perfs.seq_naive + t1 / 2;
    perfs.seq_optim = perfs.seq_optim + t2 / 2;
    perfs.par_tiled = perfs.par_tiled + t3 / 2;
    perfs.row_blocks = perfs.row_blocks + t3 / 2;
    perfs.row_blocks_bis = perfs.row_blocks_bis + t3 / 2;

    return;
}

void LCS::measure_perfs(int n){
    int i;
    for (i = 0; i < n; ++i) {
        do_perf_update();
    }

    return;
};