#include <stdlib.h>
#include <stdio.h>

#include "unique_model.h"

#include "mask_kernels.h"


void mask_density(unique_model *mod)
{
    int i;
    for (i=0;i<mod->num_ele;i++)
        mod->rho[i] = 0.0;
}

void mask_kernel(unique_model *kernel, unique_model *mask)
{
    int i;
    for (i=0;i<kernel->num_ele;i++)
    {
        kernel->epar1[i] *= mask->epar1[i];
        kernel->epar2[i] *= mask->epar1[i];
        kernel->rho[i]   *= mask->epar1[i];
    }
}

/*****************************************************************************/

void prepare_mask(double xmin, double xmax, double zmin, double zmax, unique_model *base_model, unique_model *mask)
{
    copy_unique_model(mask, base_model);

    double x, z;
    int i;
    for (i=0;i<mask->num_ele;i++)
    {
        x = mask->x[i];
        z = mask->z[i]; 

        /* default mask value is 1 */
        mask->epar1[i] = 1.0;
        mask->epar2[i] = 1.0;
        mask->rho[i]   = 1.0;

        if (x<xmin || x>xmax || z<zmin || z>zmax)
        {
            /* if outside of bounds, set mask to zero */
            mask->epar1[i] = 0.0;
            mask->epar2[i] = 0.0;
            mask->rho[i]   = 0.0;
        }
    }
}
