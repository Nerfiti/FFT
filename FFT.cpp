#include <cassert>
#include <cstdio>
#include <cstring>

#include "FFT.hpp"

static const int dec_base = 4;
static const int base = 10000;

void FFT(double _Complex *a, int size, bool inverse) 
{
    if (size&(size-1) != 0)
    {
        printf("Wrong size: %d!\n", size);
        return;
    }

    if (size == 1) {return;}

    double _Complex a0[size/2] = {};
    double _Complex a1[size/2] = {};

    for (int i = 0, j = 0; j < size/2; i+=2, ++j)
    {
        a0[j] = a[i];
        a1[j] = a[i+1];
    }

    FFT(a0, size/2, inverse);
    FFT(a1, size/2, inverse);

    double angle = 2*M_PI/size * (inverse ? -1 : 1);
    double _Complex w  = 1;
    double _Complex wn = cos(angle) + sin(angle) * I;

    for (int i = 0; i < size/2; ++i)
    {
        a[i] = (a0[i] + w * a1[i]) / (inverse ? 2 : 1);
        a[i + size/2] = (a0[i] - w * a1[i]) / (inverse ? 2 : 1);
        w *= wn;
    }
}

void Multiply_Polynoms(int *a, int first_size, int *b, int second_size, int *result)
{
    int max_size = (first_size > second_size) ? first_size : second_size;

    int size = 1;
    while (size < max_size) {size <<= 1;}
    size <<= 1;

    double _Complex a_resized[size] = {};
    double _Complex b_resized[size] = {};

    for (int i = 0; i < size; ++i)
    {
        a_resized[i] = (i <  first_size) ? a[i] : 0;
        b_resized[i] = (i < second_size) ? b[i] : 0;
    }

    FFT(a_resized, size, false);
    FFT(b_resized, size, false);
    
    for (int i = 0; i < size; ++i)
    {
        a_resized[i] *= b_resized[i];
    }

    FFT(a_resized, size, true);

    for (int i = 0; i < first_size+second_size - 1; ++i)
    {
        result[i] = (int) (creal(a_resized[i]) + 0.5);
    }
}

void GetPolynom(const char *number, int *polynom, int *deg)
{  
    assert(number != nullptr);
    int len = strlen(number);

    int num_of_first_chunk_symbols = len % dec_base;
    if (num_of_first_chunk_symbols == 0) 
    {
        num_of_first_chunk_symbols = dec_base;
    }

    *deg = (int)(len / dec_base) + ((num_of_first_chunk_symbols == dec_base) ? -1 : 0);

    int value = 0;
    for (int i = 0; i < num_of_first_chunk_symbols; ++i)
    {
        value = value*10 + *number - '0';
        number++;
    }
    polynom[*deg] =  value;

    for (int i = *deg - 1; i >= 0; --i)
    {
        value = 0;
        for (int j = 0; j < dec_base; ++j)
        {
            value = value*10 + *number - '0';
            number++;
        }
        polynom[i] = value;
    } 
}

void Multiply(int *a, int first_size, int *b, int second_size, int *result)
{
    Multiply_Polynoms(a, first_size, b, second_size, result);

    int transfer = 0;
	for (int i = 0; i < first_size + second_size; ++i) 
    {
		result[i] += transfer;
		transfer = result[i] / base;
		result[i] %= base;
	}
}