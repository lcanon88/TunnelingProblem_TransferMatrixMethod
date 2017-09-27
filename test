//Including header files
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>


//Defining simulation parameters
const double end = 2.0;
const double step_size = 1E-4;
const double E_over_V0_max = 1.0;
const double V0_star_max = 10.0;
const int num_E = 10; //the number of columns (E_over_V0)
const int num_V = 20; //the number of rows (V0_star)


//Pre-declared functions
double V(double x);
gsl_complex k_fun(double x, double E_over_V0, double V0_star); //exp(-ikz)형태를 참조
double TP(double E_over_V0, double V0_star); //tunneling probability 구하는 함수
float timedifference_msec(struct timeval t0, struct timeval t1);


//Main function
int main (void)
{
    struct timeval t0;
    struct timeval t1;
    float elapsed;
    double E_over_V0, V0_star;
    int i, j;

    gettimeofday(&t0, 0);


    gsl_matrix *E_over_V0_can = gsl_matrix_calloc(num_V, num_E);
    gsl_matrix *V0_star_can = gsl_matrix_calloc(num_V, num_E);
    gsl_matrix *Tunneling_prob = gsl_matrix_calloc(num_V, num_E);

    for(i=0; i<num_V; i++) //rows
    {
        V0_star = V0_star_max*((double)(i))/((double)(num_V-1));
        for(j=0; j<num_E; j++) //columns
        {
            E_over_V0 = E_over_V0_max*((double)(j))/((double)(num_E-1));

            gsl_matrix_set(E_over_V0_can, i, j, E_over_V0);
            gsl_matrix_set(V0_star_can, i, j, V0_star);
            gsl_matrix_set(Tunneling_prob, i, j, TP(E_over_V0, V0_star) );
        }
    }

    gsl_matrix_free(E_over_V0_can);
    gsl_matrix_free(V0_star_can);
    gsl_matrix_free(Tunneling_prob);

    gettimeofday(&t1, 0);



    elapsed = timedifference_msec(t0, t1);

    //printf("E_over_V0 : %f\n", gsl_matrix_get(E_over_V0_can , 0, 0));
    //printf("V0_star : %f\n", gsl_matrix_get(V0_star_can , 0, 0));
    //printf("투과확률 : %f\n", gsl_matrix_get(Tunneling_prob , 0, 0));
    //double test = TP(0.5, 1);
    //printf("%f\n", test);

    printf("Code executed in %f milliseconds.\n", elapsed); //ms 단위로 걸린 시간 측정

    return 0;
}


//User defined function
//exp(iwt) : time-varying
//exp(-ikz) : +z-propagating plane wave
float timedifference_msec(struct timeval t0, struct timeval t1)
{
    return (t1.tv_sec - t0.tv_sec) * 1000.0f + (t1.tv_usec - t0.tv_usec) / 1000.0f;
}

double V(double x) //Potential barrier function V(x)
{
    if((-1.0 <= x) && (x <= 1.0))
    {
        return 1.0;
    }
    return 0.0;
}

gsl_complex k_fun(double x, double E_over_V0, double V0_star)
{
    return gsl_complex_conjugate( gsl_complex_sqrt_real( V0_star*(E_over_V0-V(x)) ) );
}

double TP(double E_over_V0, double V0_star)
{
    int num_x = round(2*end/step_size)+1, k;
    double R, T;
    gsl_complex x_delta, T11, T12, T21, T22, r, t;
    gsl_vector *x = gsl_vector_calloc(num_x);

    if( (E_over_V0 == 0.0)&(V0_star == 0.0) ) //물리적으로 의미없음
    {
        return 0.0;
    }
    else if( (E_over_V0 == 0.0)&(V0_star != 0.0) ) //입자에너지가 없음
    {
        return 0.0;
    }
    else if( (E_over_V0 != 0.0)&(V0_star == 0.0) ) //포텐셜 장벽이 없음
    {
        return 1.0;
    }
    else
    {

    for(k=0; k<num_x; k++) //linspace(x, 0, 1)
    {
        gsl_vector_set(x, k, -end + (double)(k)*step_size);
    }

    gsl_matrix_complex *T_left;
    gsl_matrix_complex *T_right;
    gsl_matrix_complex *temp;

    T_left = gsl_matrix_complex_calloc(2, 2);
    T_right = gsl_matrix_complex_calloc(2, 2);
    temp = gsl_matrix_complex_calloc(2, 2);

    gsl_matrix_complex_set(T_right, 0, 0, gsl_complex_rect(1,0));
    gsl_matrix_complex_set(T_right, 0, 1, gsl_complex_rect(0,0));
    gsl_matrix_complex_set(T_right, 1, 0, gsl_complex_rect(0,0));
    gsl_matrix_complex_set(T_right, 1, 1, gsl_complex_rect(1,0));

    for(k=0; k<num_x-1; k++)
    {
        x_delta = gsl_complex_sub(gsl_complex_rect(gsl_vector_get(x, k+1), 0), gsl_complex_rect(gsl_vector_get(x, k), 0));

        T11 = gsl_complex_mul( gsl_complex_mul( gsl_complex_rect(0.5, 0), gsl_complex_add( gsl_complex_rect(1, 0), gsl_complex_div(k_fun(gsl_vector_get(x, k), E_over_V0, V0_star), k_fun(gsl_vector_get(x, k+1), E_over_V0, V0_star)) ) )
                  , gsl_complex_exp( gsl_complex_mul(gsl_complex_mul(gsl_complex_rect(0, -1), k_fun(gsl_vector_get(x, k), E_over_V0, V0_star)), x_delta) ) );

        T12 = gsl_complex_mul( gsl_complex_mul( gsl_complex_rect(0.5, 0), gsl_complex_sub( gsl_complex_rect(1, 0), gsl_complex_div(k_fun(gsl_vector_get(x, k), E_over_V0, V0_star), k_fun(gsl_vector_get(x, k+1), E_over_V0, V0_star)) ) )
                  , gsl_complex_exp( gsl_complex_mul(gsl_complex_mul(gsl_complex_rect(0, 1), k_fun(gsl_vector_get(x, k), E_over_V0, V0_star)), x_delta) ) );

        T21 = gsl_complex_mul( gsl_complex_mul( gsl_complex_rect(0.5, 0), gsl_complex_sub( gsl_complex_rect(1, 0), gsl_complex_div(k_fun(gsl_vector_get(x, k), E_over_V0, V0_star), k_fun(gsl_vector_get(x, k+1), E_over_V0, V0_star)) ) )
                  , gsl_complex_exp( gsl_complex_mul(gsl_complex_mul(gsl_complex_rect(0, -1), k_fun(gsl_vector_get(x, k), E_over_V0, V0_star)), x_delta) ) );

        T22 = gsl_complex_mul( gsl_complex_mul( gsl_complex_rect(0.5, 0), gsl_complex_add( gsl_complex_rect(1, 0), gsl_complex_div(k_fun(gsl_vector_get(x, k), E_over_V0, V0_star), k_fun(gsl_vector_get(x, k+1), E_over_V0, V0_star)) ) )
                  , gsl_complex_exp( gsl_complex_mul(gsl_complex_mul(gsl_complex_rect(0, 1), k_fun(gsl_vector_get(x, k), E_over_V0, V0_star)), x_delta) ) );


        gsl_matrix_complex_set(T_left, 0, 0, T11);
        gsl_matrix_complex_set(T_left, 0, 1, T12);
        gsl_matrix_complex_set(T_left, 1, 0, T21);
        gsl_matrix_complex_set(T_left, 1, 1, T22);

        //gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), &T_left.matrix, &T_right.matrix, gsl_complex_rect(0,0), &temp.matrix); //  temp = T_left * T_right (행렬 곱)
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), T_left, T_right, gsl_complex_rect(0,0), temp); //  temp = T_left * T_right (행렬 곱)

        gsl_matrix_complex_set(T_right, 0, 0, gsl_matrix_complex_get(temp, 0, 0));
        gsl_matrix_complex_set(T_right, 0, 1, gsl_matrix_complex_get(temp, 0, 1));
        gsl_matrix_complex_set(T_right, 1, 0, gsl_matrix_complex_get(temp, 1, 0));
        gsl_matrix_complex_set(T_right, 1, 1, gsl_matrix_complex_get(temp, 1, 1));

    }

    T11 = gsl_matrix_complex_get(temp, 0, 0);
    T12 = gsl_matrix_complex_get(temp, 0, 1);
    T21 = gsl_matrix_complex_get(temp, 1, 0);
    T22 = gsl_matrix_complex_get(temp, 1, 1);

    r = gsl_complex_mul( gsl_complex_div(gsl_complex_mul(gsl_complex_rect(-1, 0), T21), T22), gsl_complex_exp(gsl_complex_mul(gsl_complex_mul(gsl_complex_rect(0, -2), k_fun(gsl_vector_get(x, 0), E_over_V0, V0_star)), gsl_complex_rect(gsl_vector_get(x, 0), 0))) );

    t = gsl_complex_mul(gsl_complex_sub(T11, gsl_complex_div(gsl_complex_mul(T12, T21), T22)), gsl_complex_exp(gsl_complex_mul(gsl_complex_rect(0,1), gsl_complex_sub(gsl_complex_mul(k_fun(gsl_vector_get(x, num_x-1), E_over_V0, V0_star), gsl_complex_rect(gsl_vector_get(x, num_x-1), 0)), gsl_complex_mul(k_fun(gsl_vector_get(x, 0), E_over_V0, V0_star), gsl_complex_rect(gsl_vector_get(x, 0), 0))))) );

    gsl_vector_free(x);
    gsl_matrix_complex_free(T_left);
    gsl_matrix_complex_free(T_right);
    gsl_matrix_complex_free(temp);

    R = gsl_complex_abs2(r);
    T = gsl_complex_abs2(t);

    if(R+T-1.0>1E-6)
    {
        return 10.0;
    }

    return T;

    }
}
