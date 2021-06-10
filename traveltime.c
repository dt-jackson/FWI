#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include "seis_data.h"
#include "gauss_conv.h"
#include "traveltime.h"



void cross_correlation_fft(double *x, int Nx, double *y, int Ny, double *h) 
{
    /* h should be pre allocated to have a length of at least Nx+Ny-1 */

    int i;
    int cc_len = Nx + Ny - 1;

    int fft_len = next_power2(cc_len);
    int freq_len = fft_len/2 + 1;

    double *x_in = fftw_malloc(sizeof(*x_in) * fft_len);
    double *y_in = fftw_malloc(sizeof(*y_in) * fft_len);
    double *h_out = fftw_malloc(sizeof(*h_out) * fft_len);
    fftw_complex *X = fftw_malloc(sizeof(*X) * freq_len);
    fftw_complex *Y = fftw_malloc(sizeof(*Y) * freq_len);

    fftw_plan Xplan, Yplan, reverse_plan;
    Xplan = fftw_plan_dft_r2c_1d(fft_len, x_in, X, FFTW_MEASURE);
    Yplan = fftw_plan_dft_r2c_1d(fft_len, y_in, Y, FFTW_MEASURE);
    reverse_plan = fftw_plan_dft_c2r_1d(fft_len, X, h_out, FFTW_MEASURE);

    
    /* copy data */
    for (i=0;i<Nx;i++)
        x_in[i] = x[i];
    for (i=Nx;i<fft_len;i++)
        x_in[i] = 0.0;

    for (i=0;i<Ny;i++)
        y_in[i] = y[i];
    for (i=Ny;i<fft_len;i++)
        y_in[i] = 0.0;


    /* perform FFT */
    fftw_execute(Xplan);
    fftw_execute(Yplan);

    /* multiply X by conjugate of Y */
    for (i=0;i<freq_len;i++)
        X[i] = X[i] * conj(Y[i]);

    /* perform reverse FFT on result */
    fftw_execute(reverse_plan);

    /* get the data we care about */
    /* positive shifts */
    for (i=0;i<Nx;i++)
        h[i] = h_out[i];
    /* negative shifts */
    for (i=0;i<Ny-1;i++)
    {
        h[Nx + i] = h_out[fft_len-Ny+1+i];
    }


    /* clean up */
    fftw_destroy_plan(Xplan);
    fftw_destroy_plan(Yplan);
    fftw_destroy_plan(reverse_plan);//
//



    fftw_free(Y);
    fftw_free(X);
    fftw_free(h_out);
    fftw_free(y_in);
    fftw_free(x_in);

    return; 
}



double traveltime_diff(double *seis_syn, double *seis_obs, int N, double dt)
{
    int cc_len = 2*N - 1;
    double *hh = malloc(sizeof(*hh) * cc_len);

    /* perform cross correlation */
    cross_correlation_fft(seis_syn, N, seis_obs, N, hh);

    /* pick max shift */
    double max_h = hh[0];
    int max_ind = 0;
    int i;
    /* positive shift */
    for (i=1;i<N;i++)
    {
        if (hh[i] > max_h)
        {
            max_h = hh[i];
            max_ind = i;
        }
    }

    /* negative shift */
    for (i=N;i<cc_len;i++)
    {
        if (hh[i] > max_h)
        {
            max_h = hh[i];
            max_ind = -(cc_len - i);
        }
    }

    free(hh);

    /* return time shift */
    return (max_ind * dt);
}


double traveltime_misfit(seis_data *obs, seis_data *syn, double *shifts)
{
    /* returns the total cross-correlation traveltime misfit         */
    /* a matrix of the shifts is returned by address                 */
    /* shifts should be pre-allocated with num_src x num_rec entries */
    
    int i, j;
    int seis_ind;
    int shift_ind;
    double total_misfit = 0.0;
    double tmp_misfit;
    double *obs_trace, *syn_trace;

    int num_src = obs->num_src;
    int num_rec = obs->num_rec;
    int nstep   = obs->num_samples;
    double dt   = obs->dt;


    for (i=0;i<num_src;i++)
    {
        for (j=0;j<num_rec;j++)
        {
            seis_ind = (i * num_rec * nstep) + (j * nstep);
            shift_ind = (i * num_rec) + j;

            /* get pointers to current trace */
            syn_trace = syn->traces + seis_ind;
            obs_trace = obs->traces + seis_ind;

            /* calculate cross-correlation traveltime */
            tmp_misfit = traveltime_diff(syn_trace, obs_trace, nstep, dt);

            /* save shift and add squared traveltime difference */
            shifts[shift_ind] = tmp_misfit;
            total_misfit += (tmp_misfit * tmp_misfit);
        }
    } 

    /* multiply total_misfit by half - see definition of traveltime misfit */
    return (0.5 * total_misfit);
}


void traveltime_adjoint(seis_data *syn, seis_data *adjoint, double *shifts)
{
    /* note that specfem will do the time reversal */

    int i, j, k;
    int seis_ind, shift_ind;
    double *syn_trace, *adj_trace;
    

    int num_src = syn->num_src;
    int num_rec = syn->num_rec;
    int nstep   = syn->num_samples;
    double dt   = syn->dt;

    double *syn_vel = malloc(sizeof(*syn_vel) * nstep);

    double Nr, qq, tmp_misfit;


    double tol = 1e-16;

    for (i=0;i<num_src;i++)
    {
        for (j=0;j<num_rec;j++)
        {
            seis_ind = (i * num_rec * nstep) + (j * nstep);
            shift_ind = (i * num_rec) + j;

            /* get pointers to current trace */
            syn_trace = syn->traces + seis_ind;
            adj_trace = adjoint->traces + seis_ind;

            /* calculate synth velocity by finite differences */
            first_difference(syn_trace, syn_vel, nstep, dt);

            /* calculate norm factor (negative integral of syn_vel squared) */
            /* note that this could also be positive integral of syn_acc*syn_disp */
            /* see Liu 2007 (think about it in the frequency domain) */
            Nr = 0.0;
            for (k=0;k<nstep;k++)
            {
                Nr += (syn_vel[k] * syn_vel[k]);
            }
            Nr *= (-dt);

            tmp_misfit = shifts[shift_ind];

            /* divide misift by norm factor - check if Nr is 0 first */            
            if (fabs(Nr) > tol) 
            {
                qq = tmp_misfit / Nr;
            }
            else
            {
                qq = 0.0;
            }

            /* get adjoint trace by  multiplying velocity by qq */
            for (k=0;k<nstep;k++)
            {
                adj_trace[k] = qq * syn_vel[k];
            }


        }
    }

    free(syn_vel);

    return;
}




void first_difference(double *disp, double *vel, int N, double dt)
{
    int i;

    double d1 = 1.0 / (2.0 * dt);
    double d2 = 1.0 / dt;

    /* central differences for the middle */
    for (i=1;i<N-1;i++)
    {
        vel[i] = d1 * (disp[i+1] - disp[i-1]);
    }


    /* forward/backward difference at beginning and end */
    vel[0]   = d2 * (disp[1] - disp[0]);
    vel[N-1] = d2 * (disp[N-1] - disp[N-2]);

    return;
}

void second_difference(double *disp, double *acc, int N, double dt)
{
    /* not fully implemented yet - right now just calls first_difference twice, which is not efficient */
    double *vel = malloc(sizeof(*vel) * N);
    first_difference(disp, vel, N, dt);
    first_difference(vel, acc, N, dt);
    free(vel);
    return;
}