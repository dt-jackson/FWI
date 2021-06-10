#ifndef CUBIC_RECONSTRUCTION_H
#define CUBIC_RECONSTRUCTION_H



double cubic_interp(double x1, double x2, double fx1, double fx2, double fpx1, double fpx2, double x);
void cubic_recon(double *data_in, double *data_out, int N, double dt, double t0, double upper_cutoff, double lower_cutoff);

#endif