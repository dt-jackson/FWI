#include <stdlib.h>
#include <stdio.h>

#ifdef MAKE_TEST
    #include <math.h>
#endif

#include "cubic_reconstruction.h"

double cubic_interp(double x1, double x2, double fx1, double fx2, double fpx1, double fpx2, double x)
{
    /* perform cubic interpolation based on two points with known       */
    /* function values and slopes                                       */
    /* x1, x2 are points with corresponding function values fx1 and fx2 */
    /* and corresponding slopes fpx1 and fpx2                           */
    /* The function will be interpolated at x                           */

    double tx = (x - x1) / (x2 - x1);
    double a =  fpx1 * (x2 - x1) - (fx2 - fx1);
    double b = -fpx2 * (x2 - x1) + (fx2 - fx1);

    double fx;
    fx  = (1 - tx) * fx1
        + tx * fx2 
        + tx * (1 - tx) * ((1 - tx) * a + tx * b);

    return fx;
}


void cubic_recon(double *data_in, double *data_out, int N, double dt, double t0, double upper_cutoff, double lower_cutoff)
{
    /* this can be performed in-place (i.e., data_in == data_out) or out of place (data_in != data_out) */

    /* maximum number of clipped regions */
    int max_regions = 100;

    double *t = malloc(sizeof(*t) * N);
    if (!t)
    {
        fprintf(stderr,"Error: cubic_recon: memory allocation failure\n");
        exit(1);
    }

    int no_regions = 0;
    int in_region = 0;
    int cstart[max_regions];
    int cend[max_regions];

    /* determined clipped regions */
    int i, j;
    for (i=0;i<N;i++)
    {
        /* calculate t */
        t[i] = (dt * i) + t0;

        /* copy data */
        data_out[i] = data_in[i];

        /* determine start of clipped region */
        if (data_out[i] > upper_cutoff || data_out[i] < lower_cutoff)
        {
            if (!in_region)
            {
                cstart[no_regions] = i;
                in_region = 1;
            }
            continue;
        }

        /* determine end of clipped region */
        if (in_region)
        {
            cend[no_regions] = i;
            in_region = 0;
            no_regions++;
        }
            
    }

    /* perform interpolation in clipped regions */
    double m1, m2;
    double x1, x2;
    double fx1, fx2;
    int clen;
    double x, fx;
    for (i=0;i<no_regions;i++)
    {
        /* point 1 */
        x1 = t[cstart[i] - 1];

        /* point 1 function value */
        fx1 = data_out[cstart[i] - 1];

        /* points 1 slope */
        m1 = (data_out[cstart[i] - 1] - data_out[cstart[i] - 2])
           / (t[cstart[i] - 1] - t[cstart[i] - 2]);


        /* point 2 */
        x2 = t[cend[i]];

        /* point 2 function value */
        fx2 = data_out[cend[i]];

        /* point 2 slope */
        m2 = (data_out[cend[i] + 1] - data_out[cend[i]])
           / (t[cend[i] + 1] - t[cend[i]]);

        /* length of clipped area */
        clen = cend[i] - cstart[i];

        /* perform interpolation for each point in clipped area */
        for (j=0;j<clen;j++)
        {
            x = t[cstart[i] + j];
            fx = cubic_interp(x1, x2, fx1, fx2, m1, m2, x);
            data_out[cstart[i] + j] = fx;
        }

    }


    free(t);

}


#ifdef MAKE_TEST
int main()
{
    int N = 200;
    double dt = 1;
    double t0 = 0;
    double d[N];
    double r[N];
    double t;

    int i;
    for (i=0;i<N;i++)
    {
        t = i*dt + t0;
        d[i] = sin(2.0 * M_PI * t * (1.0/(double)N));
    }

    double upper_cutoff = 0.5;
    double lower_cutoff = -0.8;

    cubic_recon(d, r, N, dt, t0, upper_cutoff, lower_cutoff);

    for (i=0;i<N;i++)
    {
        t=i*dt + t0;
        printf("%f\t%f\t%f\n",t,d[i],r[i]);
    }

    return 0;
}
#endif