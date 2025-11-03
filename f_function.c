#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f_functionh.c"

// ?p????? p??ß??µat??: 1 ? 2
#define PROBLEM_ID 1
// #define PROBLEM_ID 2

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double f(double t, double y) {
#if PROBLEM_ID == 1
    (void)t; // unused
    return 2.0 * y;
#elif PROBLEM_ID == 2
    (void)y; // unused
    return 1.0 - 2.0 * M_PI * sin(2.0 * M_PI * t);
#else
    fprintf(stderr, "Error: PROBLEM_ID must be 1 or 2\n");
    exit(EXIT_FAILURE);
#endif
}

double y_exact(double t) {
#if PROBLEM_ID == 1
    return exp(2.0 * t);
#elif PROBLEM_ID == 2
    return t + cos(2.0 * M_PI * t);
#else
    fprintf(stderr, "Error: PROBLEM_ID must be 1 or 2\n");
    exit(EXIT_FAILURE);
#endif
}

