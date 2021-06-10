#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include <float.h>

#include <fftw3.h>


#include "gauss_conv.h"
#include "seis_data.h"

#include "calculate_source.h"


void calculate_source(seis_data *obs, seis_data *syn, double *syn_src, double *obs_src, int src_no)
{
    double lambda1, lambda2; /* <---Regularization parameters for Tikhonov regularization */

    int npts = obs->num_samples;
    int num_rec = obs->num_rec;

    int N = npts;       /* fft size */
    int Nfreq = N/2 +1; /* size of freq vector */

    int i, j; /* counters for loops */

    /*************************************************************************/

    /* allocate space for vector of sizes */
    int *nn = malloc(sizeof(*nn) * num_rec);

    /* allocate space for syn/obs and Fourier transforms */
    double *sd = fftw_malloc(sizeof(*sd) * N * num_rec);
    fftw_complex *SD = fftw_malloc(sizeof(*SD) * Nfreq * num_rec);

    /* allocate space for syn/obs source and Fourier transforms  */
    double *k = fftw_malloc(sizeof(*k) * N);
    fftw_complex *K = fftw_malloc(sizeof(*K) * Nfreq);

    /* Allocate space for Fourier transform of Green's function (impulse response) */ 
    fftw_complex *G = fftw_malloc(sizeof(*G) * Nfreq * num_rec);

    /* allocate some temporary working space */
    fftw_complex *temp = fftw_malloc(sizeof(*temp) * 3 * N);
    fftw_complex *num = temp + N;
    fftw_complex *den = temp + 2 * N;


    /* check memory allocations */
    if (!nn || !sd || !SD || !k || !K || !G || !temp)
    {
        fprintf(stderr,"Error: calculate_source: memory allocation failure\n");
        exit(1);
    }

    /*************************************************************************/

    /* fill array of sizes - all sizes are the same */
    for (i=0;i<num_rec;i++) nn[i] = N;

    /* make FFTW plans */
    fftw_plan plan_forward_sd, plan_forward_k, plan_backward_K;

    /* synth/obs data - forward only */
    plan_forward_sd = fftw_plan_many_dft_r2c(1,nn,num_rec,sd,NULL,1,N,SD,NULL,1,Nfreq,FFTW_MEASURE);

    /* synth/obs source - forward and backward */
    plan_forward_k = fftw_plan_dft_r2c_1d(N,k,K,FFTW_MEASURE);
    plan_backward_K = fftw_plan_dft_c2r_1d(N,K,k,FFTW_MEASURE);

    /*************************************************************************/
    /* copy synth source and data over */
    for (j=0;j<npts;j++)
    {
        k[j] = syn_src[j];
    }
    for (j=npts;j<N;j++)
    {
        /* zero padding */
        k[j] = 0.0;
    }

    for (i=0;i<num_rec;i++)
    {
        for (j=0;j<npts;j++)
        {
            sd[i*N + j] = syn->traces[src_no*npts*num_rec + i*npts + j];
        }
        for (j=npts;j<N;j++)
        {
            /* zero padding */
            sd[i*N + j] = 0.0;
        }
    }

    /* Execute FFTW plans */
    fftw_execute(plan_forward_sd);
    fftw_execute(plan_forward_k);


    /*************************************************************************/

    /* calculate denominator of first fraction */
    /* zero out num and den of second fraction */
    for (j=0;j<N;j++)
    {
        /* pre computing this reduces the number of divisions by a factor of num_rec */
        temp[j] = 1.0 / (conj(K[j]) * K[j] + lambda1);
        num[j] = 0.0;
        den[j] = 0.0;
    }

    /* solve for Green's function (freq domain) */
    for (i=0;i<num_rec;i++)
    {
        /* alternatively, call spectral_division_wl or spectral_division_tik */
        /* in these cases temp does not need to be the inverse as above */
        for (j=0;j<N;j++)
        {
            G[i*N + j] = conj(K[j]) * SD[i*N + j] * temp[j];
        }
    }

    /*************************************************************************/

    /* copy over observed data and take Fourier transform */
    for (i=0;i<num_rec;i++)
    {
        for (j=0;j<npts;j++)
        {
            sd[i*N + j] = obs->traces[src_no*npts*num_rec + i*npts + j];
        }
        for (j=npts;j<N;j++)
        {
            /* zero padding */
            sd[i*N + j] = 0.0;
        }
    }

    fftw_execute(plan_forward_sd);

    /*************************************************************************/
    /* solve for obs source in least squares sense */

    for (i=0;i<num_rec;i++)
    {
        for (j=0;j<N;j++)
        {                
            num[j] += conj(G[i*N + j]) * SD[j];
            den[j] += conj(G[i*N + j]) * G[i*N + j];
        }
    }

    /* alternatively, call spectral_division_wl or spectral_division_tik */
    for (j=0;j<N;j++)
    {
        K[j] = num[j] / (den[j] +lambda2);
    }

    /* Execute FFTW plan for backwards transform of k */
    fftw_execute(plan_backward_K);

    /* scale k */
    double sf = 1.0 / (double)N;
    for (j=0;i<N;j++)
    {
        k[j] = sf * k[j];
    }


    /* get minimum phase reconstruction. write result in output vector */
    min_phase_reconstruction(k, obs_src, N);


}


void min_phase_reconstruction(const double *signal_in, double *signal_out, const int npts)
{
    /* construct a minimum phase wavelet from the data in signal_in                 */
    /* npts is the actual length of the data without any padding                    */
    /* signal_out should be allocated already and able to hold at least npts values */

    /* this restricts us to even N */
    int N = next_power2(npts*8);

    /* allocate temporary arrays */
    fftw_complex *x1 = malloc(sizeof(*x1) * N);
    fftw_complex *x2 = malloc(sizeof(*x2) * N);
    if (!x1 || !x2)
    {
        fprintf(stderr, "Error: min_phase_reconstruction: Memory allocation failure\n");
        exit(1);
    }

    /* make fftw plans */
    fftw_plan ft_plan_for, ft_plan_bac;
    ft_plan_for = fftw_plan_dft_1d(N, x1, x2, FFTW_FORWARD, FFTW_MEASURE);
    ft_plan_bac = fftw_plan_dft_1d(N, x2, x1, FFTW_BACKWARD, FFTW_MEASURE);


    int i;
    /* copy data over */
    for (i=0;i<npts;i++)
    {
        x1[i] = signal_in[i];
    }
    /* put zero padding at the end */
    for (i=npts;i<N;i++)
    {
        x1[i] = 0.0;
    }

    /* forward fft */
    fftw_execute(ft_plan_for);

    /* take log of abs */ /* <---TODO: need to handle values with 0 magnitude */
    for (i=0;i<N;i++)
    {
        x2[i] = log(cabs(x2[i]));
    }

    /* inverse fft */
    fftw_execute(ft_plan_bac);

    /* for inverse fft, need to divide values by fft length */
    double scale_fact = 1.0 / (double)N;
    int half_ind = N/2;

    /* window signal */
    /* first element (zero frequency) is unchanged - take real part and scale */
    x1[0] = scale_fact * creal(x1[0]);

    /* take real part, scale, and window - positive freqs by 2, negative freqs by 0 */
    for (i=1;i<N/2;i++)
    {
        x1[i] = 2.0 * scale_fact * creal(x1[i]); 
        x1[i + half_ind] = 0.0;
    }
    /* element at N/2 (nyquist frequency) is unchanged - take real part and scale */
    x1[half_ind] = scale_fact * creal(x1[half_ind]);

    /* execute forward fft of windowed signal */
    fftw_execute(ft_plan_for);

    /* run signal through (complex) exponential */
    for (i=0;i<N;i++)
    {
        x2[i] = cexp(x2[i]);
    }

    /* take inverse fft */
    fftw_execute(ft_plan_bac);

    /* take real part and scale data. place result in output array*/
    for (i=0;i<npts;i++) 
    {
        signal_out[i] = scale_fact * creal(x1[i]);
    }



    /* clean up */
    fftw_free(x1);
    fftw_free(x2);
    fftw_destroy_plan(ft_plan_for);
    fftw_destroy_plan(ft_plan_bac);

}




void spectral_division_wl(fftw_complex *num, fftw_complex *den, int N, double wtr, fftw_complex *out)
{
    int i;
    double tol = 2 * DBL_EPSILON;

    /* temporary to hold magnitude of den */
    double *d_mag = malloc(sizeof(*d_mag) * N);
    if (!d_mag)
    {
        fprintf(stderr, "Error: spectral_division_wl: Memory allocation failure\n");
        exit(1);
    }

    double d_amp_max = cabs(den[0]);
    /* get magnitude of den and pick max */
    for (i=0;i<N;i++)
    {
        d_mag[i] = cabs(den[i]);
        d_amp_max = d_mag[i] > d_amp_max ? d_mag[i] : d_amp_max;
    }

    double wtr_use = wtr * d_amp_max;

    fftw_complex den_use;
    for (i=0;i<N;i++)
    {
        if (d_mag[i] < tol)
        {
            /* den magnitude equal to zero */
            /* set value to wtr_use */
            den_use = wtr_use;
        }
        else if (d_mag[i] < wtr_use)
        {
            /* den magnitude less than water level */
            /* keep phase the same, but change magnitude to wtr_use */
            den_use = wtr_use * (den[i] / d_mag[i]);
        }
        else 
        {
            /* den magnitude greater than water level */
            /* no change */
            den_use = den[i];
        }

        out[i] = num[i] / den_use; 
    }

    free(d_mag);
}



void spectral_division_tik(fftw_complex *num, fftw_complex *den, int N, double lambda, fftw_complex *out)
{
    int i;
    for (i=0;i<N;i++)
    {
        out[i] = (conj(den[i]) * num[i]) / (conj(den[i]) * den[i] + lambda);
    }
}




void fft_convolution(double *x, double *h, double *y, int npts)
{
    /* convolution of two equal length signals, x and h.    */
    /* x and h are both npts long                           */
    /* the result (y) is length 2*npts-1                    */
    /* y should already be allocated                        */


    /* fft size needs to be at least 2*npts -1 */
    int N = next_power2(2 * npts);

    int fft_freq_sz = N/2 + 1;

    double *x_in = fftw_malloc(sizeof(*x_in) * N * 2);
    fftw_complex *X = fftw_malloc(sizeof(*X) * fft_freq_sz * 2);
    if (!x_in || !X)
    {
        fprintf(stderr,"Error: fft_convolution: memory allocation failure\n");
        exit(1);
    }


    double *h_in = x_in + N;
    fftw_complex *H = X + fft_freq_sz;

    int how_many = 2; 
    int n[2] = {N,N};
    int stride = 1;
    int rank = 1;
    fftw_plan fplan, bplan;
    fplan = fftw_plan_many_dft_r2c(rank,n,how_many,x_in,NULL,stride,N,X,NULL,stride,fft_freq_sz,FFTW_MEASURE);
    bplan = fftw_plan_dft_c2r_1d(N,X,x_in,FFTW_MEASURE);

    int i;
    /* copy data */
    for (i=0;i<npts;i++)
    {
        x_in[i] = x[i];
        h_in[i] = h[i];
    }
    /* add zero padding */
    for (i=npts;i<N;i++)
    {
        x_in[i] = 0.0;
        h_in[i] = 0.0;
    }

    /* execute forward plan */
    fftw_execute(fplan);


    /* perform multiplication - keep result in X */
    for (i=0;i<fft_freq_sz;i++)
    {
        X[i] *= H[i];
    }

    /* execute backward plan */
    fftw_execute(bplan);

    /* scale and copy data to output array */
    double sf = 1.0 / (double)N;
    for (i=0;i<2*npts-1;i++)
    {
        y[i] = sf * x_in[i];
    }


    /* clean up */
    fftw_destroy_plan(fplan);
    fftw_destroy_plan(bplan);

    fftw_free(x_in);
    fftw_free(X);


}








#include <lapacke.h>

void deconvolution_system(const double *x, double *h, const double *y, const int npts, const double lambda)
{
    /* NOTE: This is really only practical when npts is fairly small        */
    /*                                                                      */
    /* the convolution problem can be formulated as a matrix multiplication */
    /* X*h = y, where the columns of X contain shifted copies of x          */
    /* Therefore, the deconvolution problem is solving this overdetermined  */
    /* system for h. To do this, we solve (X'*X + lambda*I)*h = X'*y        */
    /* where ' denotes matrix transpose, I is the identity matrix, and      */
    /* lambda is a Tikhonov regularization parameter                        */
    /* x and h should be length npts and y is length (2*npts-1)             */

    int N = next_power2(2*npts);
    int freq_sz = N/2 + 1;
    int conv_sz = 2*npts -1;


    double *X = malloc(sizeof(*X) * npts * conv_sz);
    double *xx = malloc(sizeof(*xx) * 3 * npts);
    double *yy = malloc(sizeof(*yy) * conv_sz);
    double *s = malloc(sizeof(*s) * npts);  /* s holds the singular values */


    /* copy data and zero pad */
    int i, j;    
    for (i=0;i<npts;i++)
    {
        xx[i] = 0.0;
        xx[i+npts] = x[i];
        xx[i+2*npts] = 0.0;
        yy[i] = y[i];
    }
    for (i=npts;i<conv_sz;i++)
    {
        yy[i] = y[i];
    }

    for (i=0;i<npts;i++)
    {
        for (j=0;j<conv_sz;j++)
        {
            X[i*conv_sz + j] = xx[(npts - i) + j];
        }
    }


    char trans = 'N';

    int m = conv_sz; /* number of rows in X */
    int n = npts; /* number of cols in X */

    int nrhs = 1; /* number of right hand side (1) */

    int lda = m; /* leading dim of A, m for col major layout */
    int ldb = m; /* leading dim of B */

    double rcond = lambda;
    int rank;

    int info;

    info = LAPACKE_dgelss(LAPACK_COL_MAJOR, m, n, nrhs, X, lda, yy, ldb, s, rcond, &rank);

    if (info != 0)
    {
        fprintf(stderr, "Error: deconvolution_system: Lapack failure\n");
        exit(1);
    }

    /* copy solution over to h */
    for (i=0;i<npts;i++)
    {
        h[i] = yy[i];
    }


}