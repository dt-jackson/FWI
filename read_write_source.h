#ifndef READ_WRITE_SOURCE_H
#define READ_WRITE_SOURCE_H

/* read and write source from/to a two column text file */
/* when reading, the first column is not read */

void read_source(char *fname, int N, double *src);
void write_source(char *fname, double *src, int N, double t0, double dt);

/* define a Gaussian approximation to an impulse */
void define_gaussian(double *g, int N, double t0, double dt, double tshift, double fmax);


#endif