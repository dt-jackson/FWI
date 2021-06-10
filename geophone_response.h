#ifndef GEOPHONE_RESPONSE_H
#define GEOPHONE_RESPONSE_H

#include <complex.h>
#include <fftw3.h>

typedef struct geophone
{
    double SG;          /* geophone scale factor            */
    double w0;          /* geophone resonant frequency (Hz) */
    double lambda;      /* geophone damping ratio           */
} geophone;


fftw_complex geophone_transfer(double w, geophone gp, double dt, char trans_type);
void geophone_motion2volt(double *motion_in, double *volt_out, int N, double dt, geophone gp, char motion_type);
void geophone_volt2motion(double *volt_in, double *motion_out, int N, double dt, geophone gp, char motion_type);

#endif