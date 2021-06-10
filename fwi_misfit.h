#ifndef FWI_MISFIT_H
#define FWI_MISFIT_H

#include "fwi.h"

int get_ndim(void *model);
void model2vector(void *model, double *x);
void vector2model(void *model, const double *x);

void hash_model(const double *x, const int N, unsigned char *mhash);

double fwi_misfit(const double *x, const int N, void *model);
void fwi_gradient(const double *x, double *gradient, const int N, void *model);

void write_adjoint_sources(seis_data *d, double stdev);
void write_new_model(fwi *ff, const double *x);
void read_synth(seis_data *s);
void process_synth(seis_data *s); /* <---this is a placeholder function for now */
void* fwi_misfit_allocate(long nbytes);

void enforce_par_bounds(unique_model *mod, parameters *p);


void mask_src_locs(unique_model *kernel, double src_loc_x, double src_loc_z, double dist_co, double exponent);

#endif