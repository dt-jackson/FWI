#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include <fftw3.h>

#include "seis_data.h"
#include "waveform.h"

double waveform_misfit(seis_data *resid, double data_stdev)
{
    /* this is an approximation to the integral, should be very */
    /* close as long as dt is not too large                     */

    long i;
    long tot_samp = resid->num_samples * resid->num_rec * resid->num_src;
    double mf = 0;

    for (i=0;i<tot_samp;i++)
    {
        mf += resid->traces[i] * resid->traces[i];
    }

    double inv_variance = 1.0 / (data_stdev * data_stdev);
    return 0.5 * mf * inv_variance * resid->dt;
}

void calc_residual(seis_data *obs, seis_data *syn, seis_data *resid)
{
    /* make sure the data are compatible */
    if (!check_compatible(obs,syn) || !check_compatible(syn, resid))
    {
        fprintf(stderr,"Data not compatible\n");
        exit(-100);
    }

    long i;
    long tot_samp = obs->num_samples * obs->num_rec * obs->num_src;

    for (i=0;i<tot_samp;i++)
    {
        resid->traces[i] = syn->traces[i] - obs->traces[i];
    }
}


void hilbert(const double *x_in, double *h_out, const int N)
{
    int fft_size = N;

    fftw_complex *x = fftw_malloc(sizeof(*x) * fft_size);
    fftw_complex *X = fftw_malloc(sizeof(*X) * fft_size);

    fftw_plan fplan, bplan;
    fplan = fftw_plan_dft_1d(fft_size, x, X, FFTW_FORWARD, FFTW_MEASURE);
    bplan = fftw_plan_dft_1d(fft_size, X, x, FFTW_BACKWARD, FFTW_MEASURE);

    /* copy data into array */
    int i;
    for (i=0;i<N;i++)
        x[i] = x_in[i];
    for (i=N;i<fft_size;i++)
        x[i] = 0.0;

    /* perform forward FFT */
    fftw_execute(fplan);


    /* modify to get DFT of analytic signal */
    /* first element is unchanged */
    /* positive part is doubled */
    for (i=1;i<fft_size/2;i++)
        X[i] = 2.0*X[i];

    /* negative part is set to zero */
    for (i=fft_size/2+1;i<fft_size;i++)
        X[i] = 0.0;


    /* perform inverse FFT */
    fftw_execute(bplan);

    /* take imaginary part and apply scale factor */
    /* result is placed in output array */
    double sf = 1.0/((double)fft_size);
    for (i=0;i<N;i++)
        h_out[i] = sf * cimag(x[i]);


    /* clean up */
    fftw_destroy_plan(fplan);
    fftw_destroy_plan(bplan);
    fftw_free(x);
    fftw_free(X);

}



void calc_envelope(const double *x_in, double *env_out, const int N)
{

    /* calculate hilbert transform */
    hilbert(x_in, env_out, N);


    /* analytic signal is x + i*H{x} */
    /* take norm of analytic signal to get envelope */
    int i;
    for (i=0;i<N;i++)
        env_out[i] = sqrt((env_out[i] * env_out[i]) + (x_in[i] * x_in[i]));


}

double env_diff_misfit(const seis_data *s,const seis_data *d, seis_data *resid)
{

    double misfit = 0.0;


    int num_src = s->num_src;
    int num_rec = s->num_rec;
    int nstep = s->num_samples;

    double *Es = malloc(sizeof(*Es) * nstep);
    double *Ed = malloc(sizeof(*Ed) * nstep);
    double *temp1 = malloc(sizeof(*temp1) * nstep);
    double *temp2 = malloc(sizeof(*temp2) * nstep);
    double *temp3 = malloc(sizeof(*temp3) * nstep);




    int samp_per_src = num_rec*nstep;

    int i, j, k;
    double *trace_ptr_s, *trace_ptr_d, *trace_ptr_r;
    double Erat;

    double tol = 1e-30;

    for (i=0;i<num_src;i++)
    {
        for (j=0;j<num_rec;j++)
        {
            /* get trace pointers */
            trace_ptr_s = s->traces     + samp_per_src*i + nstep*j;
            trace_ptr_d = d->traces     + samp_per_src*i + nstep*j;
            trace_ptr_r = resid->traces + samp_per_src*i + nstep*j;

            if (zero_trace(trace_ptr_s, nstep, tol))
            {
                /* set adjoint souce to zero */
                for (k=0;k<nstep;k++)
                {
                    trace_ptr_r[k] = 0.0;
                }

                /* move on to next trace */
                continue;
            }


            /* calculate envelopes */
            calc_envelope(trace_ptr_s, Es, nstep);
            calc_envelope(trace_ptr_d, Ed, nstep);

            /* calculate hilbert transform of synth for adjoint source */ 
            hilbert(trace_ptr_s, temp1, nstep);

            /* sum envelope difference squared  for misfit */
            for (k=0;k<nstep;k++)
            {
                misfit += ((Es[k] - Ed[k]) * (Es[k] - Ed[k]));

                /* need to make sure Es is not zero before division */
                Erat = (Es[k] - Ed[k]) / Es[k];

                /* first term of adjoint source */
                temp3[k] = Erat * trace_ptr_s[k];

                /* inner part of second term of adjoint source */
                temp2[k] = Erat * temp1[k];
            }

            /* second term of adjoint source */
            hilbert(temp2, temp1, nstep); 

            /* calculate adjoint source */
            for (k=0;k<nstep;k++)
            {
                trace_ptr_r[k] = temp3[k] - temp1[k];
            }

        }
    }

    free(temp1);
    free(temp2);
    free(temp3);
    free(Es);
    free(Ed);



    return (0.5 * misfit);


}

int zero_trace(const double *trace, const int N, const double tol)
{
    /* returns 1 if trace is completely zeroed */
    /* returns 0 otherwise */

    int i;
    for (i=0;i<N;i++)
    {
        /* as soon as we hit a non-zero value, return */
        if (fabs(trace[i]) > tol)
            return 0;
    }

    return 1;

}




int check_compatible(seis_data *obs, seis_data *syn)
{
    if (!obs || !syn)
    {
        fprintf(stderr,"Error: seis_data, check_comparible: "
                       "NULL seismic data\n");
        return 0;
    }

    int value = 1;
    /* check number of sources */
    if (obs->num_src != syn->num_src)
    {
        fprintf(stderr,"Error: seis_data, check_compatible: "
                       "Number of sources does not match\n");
        value = 0;
    }

    /* check number of receivers */
    if (obs->num_rec != syn->num_rec)
    {
        fprintf(stderr,"Error: seis_data, check_compatible: "
                       "Number of receivers does not match\n");
        value = 0;
    }

    /* check number of samples */
    if (obs->num_samples != syn->num_samples)
    {
        fprintf(stderr,"Error: seis_data, check_compatible: "
                       "Number of samples does not match\n");
        value = 0;
    }

    /* check time step */
    if (obs->dt != syn->dt)
    {
        fprintf(stderr,"Error: seis_data, check_compatible: "
                       "Time step does not match\n");
        value = 0;
    }

    return value;
}