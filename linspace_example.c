//Including header files
#include <stdio.h>
#include <stdlib.h>


double *linspace(double a, double b, int n);


int main(void)
{
    int i;
    int j = 1;

    double *test = linspace(2, 4, j);

    for(i=0; i<j; i++)
    {
        printf("%.18f\n", test[i]);
    }

    return 0;
}


double *linspace(double a, double b, int n)
{
    if(n==1)
    {
        double *x = malloc(sizeof(double) *n);
        x[0]  = a;
        return x;
    }

    int k;
    double *x = malloc(sizeof(double) *n);
    double del = (b-a) / ((double)(n-1));

    for(k=0; k<n; k++)
    {
        x[k] = a + del*(double)(k);
    }

    return x;
}
