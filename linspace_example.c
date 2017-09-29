//Including header files
#include <stdio.h>
#include <stdlib.h>


double *linspace(double a, double b, int n);


int main(void)
{
    //double *test = linspace(0.0, 1.0, 10);
    int i;

    for(i=0; i<10; i++)
    {
        printf("%.18f\n", linspace(0.0, 1.0, 10)[i]);
    }

    return 0;
}


double *linspace(double a, double b, int n)
{
    int k;
    double *x = malloc(sizeof(double) *n);
    double del = (b-a) / ((double)(n-1));

    for(k=0; k<n; k++)
    {
        x[k] = a + del*(double)(k);
    }

    return x;
}
