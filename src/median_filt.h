#ifndef MEDIAN_FILTER_H
#define MEDIAN_FILTER_H

void median_filter1d(double *x, int N, int filt_size);
void median_filter2d(double *x, int nx, int ny, int fsize_x, int fsize_y);
void median_filter3d(double *x, int nx, int ny, int nz, int fsize_x, int fsize_y, int fsize_z);

#endif