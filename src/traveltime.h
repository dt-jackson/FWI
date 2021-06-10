#ifndef TRAVELTIME_H
#define TRAVELTIME_H

#include "seis_data.h"

/* function to perform cross-correlation by FFT */
void cross_correlation_fft(double *x, int Nx, double *y, int Ny, double *h);

/* traveltime functions */
double traveltime_diff(double *seis_syn, double *seis_obs, int N, double dt);
double traveltime_misfit(seis_data *obs, seis_data *syn, double *shifts);
void traveltime_adjoint(seis_data *syn, seis_data *adjoint, double *shifts);

/* helper functions to calculate velocity and acceleration from displacement */
void first_difference(double *disp, double *vel, int N, double dt);
void second_difference(double *disp, double *acc, int N, double dt);




#endif