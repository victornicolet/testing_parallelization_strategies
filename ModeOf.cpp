//
// Created by nicolet on 26/01/18.
//

#include "ModeOf.h"

_mode_ _getmode(int* A, long n){
    _mode_ aux[n];
    _mode_ res;
    long i;

    for(i = 0 ; i< n; i++){
        aux[i] = {a[i],1};
    }


    delete [] aux;
    return res;
}

_mode_ _getmode_sorted_input(int* A, long n){
    int mode;
    long cur_cnt;
    long cnt = 0;
    long i;
    int prev = a[0];

    index = 0;
    for(i = 0; i<n; i++){
        // If we see a new element in the sorted array...
        if(prev != a[i]){
            // the old one might be the mode, so ...
            if(counts[index] > cnt){
                cnt = cur_cnt;
                mode = prev;
            }
            // .. and then we can start counting for this element
            cur_cnt = 0;
        }
        cur_cnt++;
    }
    delete [] counts;
    return {mode, cnt};
}

_mode_ _getmode_lookahead(int *A, long n){
    int mode;
    long freq;
    int B[n] = {0};
    long i = 0;
    long next = 1;

    for(i = 0; i< n; i++){
        if(B[i] = 0){ // Element not counted
            cnt = 0;
            for(long j = i; j < n; j++){
                if(A[i] == A[j]) {
                    cnt++;
                    B[j] = 1; // mark as counted
                }
            }
            if(cnt > freq) {
                freq = cnt;
                mode = A[i];
            }
        }
    }

    return {mode, freq};
}


_mode_ ModeOf::getmode() {
    return _getmode(A,n);
}