#ifndef CHEN_CLOSURE_CORE_H
#define CHEN_CLOSURE_CORE_H
extern "C" int closure_rows(const double* input, const unsigned char* force,
                           unsigned char* output, int rows, int m,
                           double alpha, int mode);
#endif
