#ifndef WAVEFORM_H
#define WAVEFORM_H

#include "seis_data.h"

int check_compatible(seis_data *obs, seis_data *syn);

/* calculate waveform residual */
void calc_residual(seis_data *obs, seis_data *syn, seis_data *resid);

/* calculate waveform misfit */
double waveform_misfit(seis_data *resid, double data_stdev);


/* calculate envelope difference misfit and adjoint source */
double env_diff_misfit(const seis_data *s,const seis_data *d, seis_data *resid);



/* calculate hilbert transform of a signal */
void hilbert(const double *x_in, double *h_out, const int N);

/* calculate envelope of a signal (magnitude of analytic signal) */
void calc_envelope(const double *x_in, double *env_out, const int N);

/* determine if a trace is completely zeroed. returns 1 if so, 0 otherwise */
int zero_trace(const double *trace, const int N, const double tol);



#endif