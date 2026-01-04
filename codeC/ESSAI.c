#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <R.h>
//#include <Rmath.h>

void coucou(int nb)
{
    int i;
    for(i = 0; i <nb; i++)
    printf("coucou", i);
}


void Affichage(int *tab, int nb)
{
    int i;
    for(i = 0; i <nb; i++)
    printf("tab[%d] = %d\n", i, tab[i]);
}


void is_dominant(double *x, double *y, int *ndim, int *set, int *length, int * result)
{
    int i, j;
    int s = 0;

printf("longueur = %d\n",*length);
printf("y = (%f,%f)\n",*(y),*(y+1));
printf("dimension = %d\n",*ndim);

    if ( set == 2)
    {
        printf("set = fail\n");
    }
    else if ( set == 1)
    {
        printf("set = safe\n");
    }
    else
    {
        printf("ERROR : set must to be safe or fail.");
        return;
    }

    for(i = 0; i < *length; i++)
    {
        s = 0;

        for(j = 0; j < *ndim; j++)
        {


           if( set == 2)
            {
                if( *(x+i + j*(*length)) >= *(y + j) )
                {
                  s += 1;
                }
            }

            else
            {
                if( *(x + i + j*(*length)) <= *(y + j))
                {
                    s += 1;
                }
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
         printf("result[%d] = %d,  x[%d] =  (%f,%f)\n", i, result[i], i+1 , *(x+i),*(x+i+(*length)) );
    }
    result = result;

}

int main(void)
{

  double x[28] = {0.4,0.1,0.5,0.6,0.6,0.8,0.8,0.9,0.9,0.4,0.5,0.2,0.3,0.2,0.3,0.4,0.5,0.5,0.4,0.2,0.4,0.3,0.1,0.5,0.6,0.6,0.7,0.9};
  double y[2] = {0.5,0.5};
  int ndim = 2;
  int set = 2;
  int length;
  length = sizeof(x)/sizeof(double)/ndim;
  int result[length];
  is_dominant( x, y, &ndim, set, &length, &result);
  return 0;
}


