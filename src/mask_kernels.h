#ifndef MASK_KERNELS_H
#define MASK_KERNELS_H

#include "unique_model.h"

typedef struct mask_parameters
{
    int no_regions;

} mask_par;



void mask_density(unique_model *mod);
void mask_kernel(unique_model *kernel, unique_model *mask);


void prepare_mask(double xmin, double xmax, double zmin, double zmax, unique_model *base_model, unique_model *mask);

#endif