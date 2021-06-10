#ifndef SEIS_DATA_H
#define SEIS_DATA_H

#include <stdio.h>

typedef struct seis_data
{
    int num_src;      /* number of sources */
    int num_rec;      /* number of receiverse per source */
    int num_samples;  /* number of samples in each trace */
    double *src_loc;   /* source location */
    double *rec_loc;  /* receiver location - array length num_traces */
    double *traces;   /* seismic data array length num_traces * num_samples */
    double dt;        /* time step */
    double delay;     /* delay/t0 (seconds) */

} seis_data;


void allocate_seis_data(seis_data *t, int num_src, int num_rec, int num_samp);
void free_seis_data(seis_data *t);
void copy_seis_data(seis_data *dest, seis_data *src);
void print_seis(FILE *ofile, seis_data *d);

#endif