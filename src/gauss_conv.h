#ifndef GAUSS_CONV_H
#define GAUSS_CONV_H


/* find the next highest power of two for 32-bit integers */
unsigned int next_power2(unsigned int v);

/* perform 1d convolution with a gaussian with variance t */
void gauss_conv1d(const double *in,double *out,int N,double t);


/* These functions are set up for row major arrays with ndim1 columns */
/* and ndim2 rows (i.e. the dim with contiguous memory goes first and */
/* is called ndim1)                                                   */
/* since gauss_conv2d is symmetric, you can just switch ndim1 and     */
/* ndim2 if you have column major                                     */

void transpose_array(double *in, double *in_transp,int ndim1, int ndim2);
void gauss_conv_array(const double *in, double *out, int n1, int n1_ele, double t1);
void gauss_conv2d(double *in, double *out, int ndim1, int ndim2, double t1, double t2);

#endif