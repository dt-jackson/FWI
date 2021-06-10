#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include <fftw3.h>

#include "geophone_response.h"


fftw_complex geophone_transfer(double w, geophone gp, double dt, char trans_type)
{
    double Fnyq = 1.0 / (2.0 * dt);     /* nyquist freq         */
    double w0 = M_PI * gp.w0 / Fnyq;    /* resonant frequency   */

    fftw_complex num, den;

    /* denomenator is the same for all three transfer functions */
    den = -(w * w) + (2 * csqrt(-1.0) * gp.lambda * w0 * w) + (w0 * w0);

    
    if (trans_type == 'd' || trans_type == 'D')
    {
        /* displacement */
        num = -gp.SG * csqrt(-1.0) * w * w * w;
    }
    else if (trans_type == 'v' || trans_type == 'V')
    {
        /* velocity */
        num = -gp.SG * w * w;
    }
    else if (trans_type == 'a' || trans_type == 'A')
    {
        /* acceleration */
        num = gp.SG * csqrt(-1.0) * w;
    }
    else
    {
        /* error */
        fprintf(stderr, "Error: geophone_transfer: transfer type '%c' unknown",trans_type);
        exit(1);
    }


    return num / den;
}





void geophone_motion2volt(double *motion_in, double *volt_out, int N, double dt, geophone gp, char motion_type)
{

    int fft_sz = N;
    int freq_sz = fft_sz / 2 + 1;

    /* allocate space for ffts */
    double *x = fftw_malloc(sizeof(*x) * fft_sz);
    fftw_complex *X = fftw_malloc(sizeof(*X) * freq_sz);


    /* make fft plans */
    fftw_plan fplan, bplan;
    fplan = fftw_plan_dft_r2c_1d(fft_sz, x, X, FFTW_MEASURE);
    bplan = fftw_plan_dft_c2r_1d(fft_sz, X, x, FFTW_MEASURE);

    /* copy data */
    int i;
    for (i=0;i<N;i++)
        x[i] = motion_in[i];

    /* zero padding */
    for (i=N;i<fft_sz;i++)
        x[i] = 0.0;

    /* execute forward fft */
    fftw_execute(fplan); 

    
    fftw_complex H;
    double w;
    double df = M_PI / (freq_sz - 1);    

    for (i=0;i<freq_sz;i++)
    {
        /* calculate frequency */
        w = df * i; 

        /* calculate transfer function */
        H = geophone_transfer(w, gp, dt, motion_type);

        /* multiply data by transfer */
        X[i] *= H;
    }

    /* inverse fft */
    fftw_execute(bplan);

    /* scale data and copy it back over */
    double sf = 1.0 / (double)fft_sz;
    for (i=0;i<N;i++)
        volt_out[i] = sf*x[i];

    /* clean up */
    fftw_destroy_plan(fplan);
    fftw_destroy_plan(bplan);
    fftw_free(X);
    fftw_free(x);

}


void geophone_volt2motion(double *volt_in, double *motion_out, int N, double dt, geophone gp, char motion_type)
{
    /* set fft sizes */
    int fft_sz = N;
    int freq_sz = fft_sz / 2 + 1;

    /* allocate memory for ffts and transfer function */
    double *x = fftw_malloc(sizeof(*x) * fft_sz);
    fftw_complex *X = fftw_malloc(sizeof(*X) * freq_sz);
    fftw_complex *H = fftw_malloc(sizeof(*H) * freq_sz);


    /* set up fft plans */
    fftw_plan fplan, bplan;
    fplan = fftw_plan_dft_r2c_1d(fft_sz, x, X, FFTW_MEASURE);
    bplan = fftw_plan_dft_c2r_1d(fft_sz, X, x, FFTW_MEASURE);

    /* copy data */
    int i;
    for (i=0;i<N;i++)
        x[i] = volt_in[i];

    /* zero padding */
    for (i=N;i<fft_sz;i++)
        x[i] = 0.0;

    /* forward fft */
    fftw_execute(fplan); 

    double w;
    double df = M_PI / (freq_sz - 1);    
    double max_H = -1;
    for (i=0;i<freq_sz;i++)
    {
        /* calculate frequency */
        w = df * i; 

        /* calculate transfer function */
        H[i] = geophone_transfer(w, gp, dt, motion_type);

        /* pick max magnitude of transfer function */
        max_H = (cabs(H[i]) > max_H) ? cabs(H[i]) : max_H;
    }

    double wl = 0.01 * max_H; /* factor of 0.01 can be adjusted if needed */
    double mag_H;
    fftw_complex phase_H, HH;

    for (i=0;i<freq_sz;i++)
    {
        mag_H = cabs(H[i]);
        if (mag_H == 0)
        {
            /* magnitude is zero, set transfer function to wl */
            HH = wl;
        }
        else if (mag_H < wl)
        {
            /* magnitude is less than wl. keep phase, but set magnitude to wl */
            phase_H = H[i] / mag_H;
            HH = wl * phase_H;
        }
        else
        {
            /* magnitude greater than wl. do nothing */
            HH = H[i];
        }

        /* spectral division */
        X[i] /= HH;
    }


    /* inverse fft */
    fftw_execute(bplan);

    /* copy data over */
    double sf = 1.0 / (double)fft_sz;
    for (i=0;i<N;i++)
        motion_out[i] = sf * x[i];

    /* clean up */
    fftw_destroy_plan(fplan);
    fftw_destroy_plan(bplan);
    fftw_free(x);
    fftw_free(X);
    fftw_free(H);


}
