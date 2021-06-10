#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include <fftw3.h>

#include "gauss_conv.h"




/*****************************************************************************/

void gauss_conv_array(const double *in, double *out, int n1, int n1_ele, double t1)
{
    /* performs a series of 1d convolutions with Gaussian  */
    /* to filter other dim, transpose array and call again */
    /* n1 is number of rows/cols, n1_ele is number of      */
    /* elements per row/col. t1 is variance of Gaussian    */ 
    int i,j;

    /* size of fft including zero padding */
    //int fft1_size = (int)next_power2((unsigned int)n1_ele);
    //fft1_size = (int)next_power2((unsigned int)fft1_size+1);
    int fft1_size = 2 * n1_ele;


    /* size of fft result from 0 to pi. Because input signal is real, */
    /* that is all we need. Symmetry can be used to get the rest      */
    int fft_freq_size = fft1_size/2 + 1;

    /* allocate local variables */
    double *in_copy = fftw_malloc(sizeof(*in_copy) * fft1_size * n1);
    fftw_complex *in_ft = fftw_malloc(sizeof(*in_ft) * fft_freq_size * n1);
    double *out_copy = fftw_malloc(sizeof(*out_copy) * fft1_size * n1);

    /* make sure memory allocation was okay */
    if (in_copy == NULL || out_copy == NULL || in_ft == NULL)
    {
        fprintf(stderr,"Memory allocation failure in gauss_fft1d\n");
        exit(1);
    }


    /* make fftw plans - since FFTW_MEASURE is used, this will take a little */
    /* longer the first time it is run, but subsequent plans should be       */
    /* quicker due to the "wisdom" of fftw                                   */
    fftw_plan ft_plan_forward;
    fftw_plan ft_plan_reverse;
    ft_plan_forward = fftw_plan_many_dft_r2c(1,&fft1_size,n1,in_copy,NULL,1,fft1_size,in_ft,NULL,1,fft_freq_size,FFTW_MEASURE);
    ft_plan_reverse = fftw_plan_many_dft_c2r(1,&fft1_size,n1,in_ft,NULL,1,fft_freq_size,out_copy,NULL,1,fft1_size,FFTW_MEASURE);


    /* initialize arrays */
    for (j=0;j<n1;j++)
    {
        for (i=0;i<n1_ele;i++)
        {
            /* copy array in */
            in_copy[i+j*fft1_size] = in[i+j*n1_ele];

            /* symmetric padding */
            //in_copy[i+j*fft1_size + n1_ele] = in[j*(n1_ele+1) - 1 - i];

            out_copy[i+j*fft1_size] = 0.0;
            //out_copy[i+j*fft1_size + n1_ele] = 0.0;
            
        }
        int cc = 0;
        for (i=n1_ele;i<fft1_size;i++)
        {

            /* use this line for zero padding */ 
            //in_copy[i+j*fft1_size] = 0;


            /* symmetric padding at the end */
            in_copy[i+j*fft1_size] = in[j*n1_ele-1-cc];
            out_copy[i+j*fft1_size] = 0;
            cc++;
        }
    }

    /* execute forward fft */
    fftw_execute(ft_plan_forward);

    /* calculate dtft values of discrete gaussian kernel */
    /* and multiply by signal                            */
    double *T = malloc(sizeof(*T)*fft_freq_size);
    double freq;
    double df = M_PI / (double)(fft_freq_size-1);
    for (i=0;i<fft_freq_size;i++)
    {
        freq = (double)i * df;
        T[i] = exp(t1 * (cos(freq)-1));
    }

    for (j=0;j<n1;j++)
    {
        for (i=0;i<fft_freq_size;i++)
        {
            in_ft[i+j*fft_freq_size] *= T[i];
        }
    }

    /* execute reverse fft */
    fftw_execute(ft_plan_reverse);     



    /* take first n1_ele values for each row and scale by fft1_size */
    double scale_factor = 1.0 / fft1_size;
    for (j=0;j<n1;j++)
    {
        for (i=0;i<n1_ele;i++)
        {
            out[i+j*n1_ele] = scale_factor * out_copy[i+j*fft1_size];
        }
    }

    
    /* free local variables */
    free(T);
    fftw_free(out_copy);
    fftw_free(in_ft);
    fftw_free(in_copy);
}

/*****************************************************************************/

void gauss_conv1d(const double *in,double *out,int N,double t)
{
    /* calculates the convolution of a 1d signal with a discrete gaussian    */
    /* kernel with variance t. Convolution is implemented in the frequency   */
    /* domain using fftw                                                     */


    int i;

    /* size of fft (including zero padding) */
    //int fft_size = (int)next_power2((unsigned int)N);
    //fft_size = (int)next_power2((unsigned int)fft_size+1);
    int fft_size = N*2;
    
    /* size of fft result from 0 to pi. Because input signal is real, */
    /* that is all we need. Symmetry can be used to get the rest      */
    int fft_freq_size = fft_size/2 + 1;

    /* allocate local variables */ 
    double *in_copy = fftw_malloc(sizeof(*in_copy)*fft_size);
    fftw_complex *in_ft = fftw_malloc(sizeof(*in_ft)*fft_freq_size);
    double *out_copy = fftw_malloc(sizeof(*out_copy)*fft_size);

    /* make sure memory allocation was okay */
    if (in_copy == NULL || out_copy == NULL || in_ft == NULL)
    {
        fprintf(stderr,"Memory allocation failure in gauss_fft1d\n");
        exit(1);
    }


    /* make fftw plans - since FFTW_MEASURE is used, this will take a little */
    /* longer the first time it is run, but subsequent plans should be       */
    /* quicker due to the "wisdom" of fftw                                   */
    fftw_plan ft_plan_forward;
    fftw_plan ft_plan_reverse;
    ft_plan_forward = fftw_plan_dft_r2c_1d(fft_size,in_copy,in_ft,FFTW_MEASURE);
    ft_plan_reverse = fftw_plan_dft_c2r_1d(fft_size,in_ft,out_copy,FFTW_MEASURE);

    /* initialize arrays */
    for (i=0;i<N;i++)
    {
        in_copy[i] = in[i];
        out_copy[i] = 0;
    }

    /* symmetric padding at the end */
    int cc = 0;
    for (i=N;i<fft_size;i++)
    {
        //in_copy[i] = 0;
        in_copy[i] = in[N-1-cc];
        out_copy[i] = 0;
        cc++;
    }

    /* execute forward fft */
    fftw_execute(ft_plan_forward);

    /* calculate dtft values of discrete gaussian kernel */
    double T;
    double freq;
    double df = M_PI / (double)(fft_freq_size - 1);
    for (i=0;i<fft_freq_size;i++)
    {
        freq = (double)i * df;
        T = exp(t * (cos(freq)-1));
        in_ft[i] *= T;
    }

    /* execute inverse fft */
    fftw_execute(ft_plan_reverse);
    
    /* Take first N values and scale by fft_size */
    double scale_factor = (1.0/fft_size);
    for (i=0;i<N;i++)
    {
        out[i] = scale_factor * out_copy[i]; 
    }

    /* free local variables */
    fftw_destroy_plan(ft_plan_reverse);
    fftw_destroy_plan(ft_plan_forward);
    fftw_free(in_ft);
    fftw_free(in_copy);
    fftw_free(out_copy);
}

/*****************************************************************************/


void gauss_conv2d(double *in, double *out, int ndim1, int ndim2, double t1, double t2)
{
    /* allocate a temporary result */
    double *temp = malloc(sizeof(*temp)*ndim1*ndim2);
    double *temp_transp = malloc(sizeof(*temp_transp)*ndim1*ndim2);
    if (!temp || !temp_transp)
    {
        fprintf(stderr, "Error: gauss_conv2d: memory allocation failure");
        exit(1);
    }

    /*************************************************************************/
    /* first method - this one works, but is a little more expensive         */

    /* smooth first dim */
    int i;
    for (i=0;i<ndim2;i++)
    {
        gauss_conv1d(in + i*ndim1, temp + i*ndim1, ndim1, t1);
    }
    
    /* transpose */
    transpose_array(temp, temp_transp, ndim1, ndim2);

    /* smooth second dim */
    for (i=0;i<ndim1;i++)
    {
        gauss_conv1d(temp_transp + i*ndim2, temp + i*ndim2, ndim2, t2);
    }

    /* transpose */
    transpose_array(temp, out, ndim2, ndim1);
    
    /*************************************************************************/
    /* second method - there is an error somewhere in gauss_conv_array,      */
    /* I suspect somewhere in the symmetric padding, but maybe in the        */
    /* fftw plan                                                             */

    ////* perform convolution along 1st dimension */
    ///gauss_conv_array(in,temp,ndim2,ndim1,t1);

    ////* transpose the result */
    ///transpose_array(temp,temp_transp,ndim1,ndim2);

    ////* perform convolution along 2nd dimension */
    ///gauss_conv_array(temp_transp,temp,ndim1,ndim2,t2);

    ////* transpose back to original */
    ///transpose_array(temp,out,ndim2,ndim1);

    /*************************************************************************/

    free(temp_transp);
    free(temp);
}

/*****************************************************************************/


unsigned int next_power2(unsigned int v)
{
    /* returns next highest power of 2 for a 32 bit integer                  */
    /* see                                                                   */
    /* https://web.archive.org/web/20160703165415/ \                         */
    /* https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2  */
    /* (should all be one url)                                               */
    /* For positive 32-bit integers. Note that an input of 0 will            */
    /* return 0 which is not a power of 2                                    */

    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return ++v;
}


/*****************************************************************************/


void transpose_array(double *in, double *in_transp,int ndim1, int ndim2)
{
    int i,j;

    for (i=0;i<ndim1;i++)
    {
        for (j=0;j<ndim2;j++)
        {
            in_transp[j+i*ndim2] = in[i+j*ndim1];
        }
    }
}