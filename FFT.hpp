#ifndef FFT_HPP
#define FFT_HPP

#include <complex.h>

void FFT (double _Complex *a, int size, bool reverse);

void Multiply_Polynoms(int *a, int first_size, int *b, int second_size, int *result);

void GetPolynom(const char *number, int *polynom, int *deg);

void Multiply(int *a, int first_size, int *b, int second_size, int *result);

#endif //FFT_HPP