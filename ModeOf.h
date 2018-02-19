//
// Created by nicolet on 26/01/18.
//

#ifndef PARALLEL_STRATEGIES_TESTING_MODEOF_H
#define PARALLEL_STRATEGIES_TESTING_MODEOF_H

struct _mode_ {
    int mode;
    long freq;
};

class ModeOf {
    int* A;
    int mode;
    long freq;
    long n;

public:
    ModeOf(int* _A, long _n): A(_A), n(_n), mode(0), freq(0) {}
    _mode_ getmode();
};


#endif //PARALLEL_STRATEGIES_TESTING_MODEOF_H
