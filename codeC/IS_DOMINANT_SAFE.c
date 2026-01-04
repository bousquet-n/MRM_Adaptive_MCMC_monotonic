#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>

void is_dominant_safe(double *x, double *y, int *ndim, int *length, int * result)
{
    int i, j;
    int s = 0;

    for(i = 0; i < *length; i++)
    {
        s = 0;
        for(j = 0; j < *ndim; j++)
        {
            if( *(x + i + j*(*length)) <= *(y + j) )
            {
                s += 1;
            }
        }
        if(s == *ndim)
        {
            result[i] = 1;
        }
        else
        {
            result[i] = 0;
        }
    }
}


