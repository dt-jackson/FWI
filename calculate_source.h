#ifndef CALCULATE_SOURCE_H
#define CALCULATE_SOURCE_H

#include <fftw3.h>






/* spectral division with water level regularization */
void spectral_division_wl(fftw_complex *num, fftw_complex *den, int N, double wtr, fftw_complex *out);

/* spectral division with tikhonov regularization */
void spectral_division_tik(fftw_complex *num, fftw_complex *den, int N, double lambda, fftw_complex *out);

/* create a minimum phase signal */
void min_phase_reconstruction(const double *in, double *out, const int npts);

/* convolution of two equal length signals */
void fft_convolution(double *x, double *h, double *y, int npts);



void deconvolution_system(const double *x, double *h, const double *y, const int npts, const double lambda);

#endif