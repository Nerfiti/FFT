#include <cstdio>

#include "FFT.hpp"

int main()
{
    const int first_size = 10;
    const int second_size = 10;
    const int result_size = first_size + second_size - 1;

    int a[first_size] = {};
    int b[second_size] = {};
    int c[result_size] = {};

    for (int i = 0; i < first_size; ++i)
    {
        a[i] = i;
    }
    for (int i = 0; i < second_size; ++i)
    {
        b[i] = 1;
    }

    Multiply_Polynoms(a, first_size, b, second_size, c);

    printf("C(x) = %d ", c[0]);
    for (int i = 1; i < result_size; ++i)
    {
        printf("+ %d*(x^%d)", c[i], i);
    }
    printf("\n");
}