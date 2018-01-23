#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

//
// Created by nicolet on 19/01/18.
//


int levenshtein(char* X, char *Y, m){
    int **D;
    int i,j;

    D = malloc(sizeof(D)*(n + 1));
    for(i = 0; i< n + 1; i++){
        D[i] = malloc(sizeof(D[i])*(n + 1));
    }

    for(i = 0; i < n+1; i++){
        D[0][i] = i;
        D[i][0] = i;
    }
/* Copyright (C) 1991-2016 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */
/* wchar_t uses Unicode 8.0.0.  Version 8.0 of the Unicode Standard is
   synchronized with ISO/IEC 10646:2014, plus Amendment 1 (published
   2015-05-15).  */
/* We do not support C11 <threads.h>.  */
      int t1, t2, t3;
     int lb, ub, lbp, ubp, lb2, ub2;
     register int lbv, ubv;
    /* Start of CLooG code */
    if (n >= 1) {
      for (t1=1;t1<=n;t1++) {
        for (t2=1;t2<=n;t2++) {
          offs = X[t1-1] == Y[t2-1] ? 2 : 0;;
          D[t1][t2] = D[t1-1][t2] + D[t1][t2-1] + D[t1-1][t2-1] + offs;;
        }
      }
    }
/* End of CLooG code */
    int res = D[n][n];

    for(i = 0; i< n; i++){
        free(D[i]);
    }
    free(D);

    return res;
}

