#ifndef FWI_H
#define FWI_H

#include "seis_data.h"
#include "ibool_model.h"
#include "unique_model.h"
#include "bilinear_interp.h"

#include "read_parameters.h"





typedef struct fwi
{
    parameters p_;
    parameters *p;

    /* observed data */
    seis_data obs_;
    seis_data *obs;

    /* synthetic data */
    seis_data syn_;
    seis_data *syn;

    /* residuals */
    seis_data res_;
    seis_data *res;

    /* ibool model - vp/vs */
    ibool_model imod_;
    ibool_model *imod;

    /* ibool conversion struct */
    convert_ibool_unique i2u_;
    convert_ibool_unique *i2u;

    /* unique model */
    unique_model umod_;
    unique_model *umod;

    /* unique model - reg grid */
    unique_model umod_reg_;
    unique_model *umod_reg;

    /* bilinear interp struct - sim2reg */
    interp2d_t *sim2reg;

    /* unique_model - reg_grid filt */
    unique_model umod_filt_;
    unique_model *umod_filt;

    /* unique model - inv grid */
    unique_model umod_inv_;
    unique_model *umod_inv;

    /* bilinear interp struct - inv2sim */
    interp2d_t *inv2sim;

    /* unique model - inv grid (kernels) */
    unique_model kernel_inv_;
    unique_model *kernel_inv;


    unique_model mask_;
    unique_model *mask;

} fwi;







void initialize_fwi(fwi *ff, parameters *p);
void free_fwi(fwi *ff);

/* functions for determining grid sizes */
int get_dec_factor(parameters *p);
void get_grid_size(parameters *p,int *nx_reg, int *nz_reg, int *nx_inv, int *nz_inv);
void make_reg_grid(unique_model *reg, double dx, double xmin, double xmax, double zmin, double zmax);
void decimate_grid(unique_model *inv, unique_model *reg, int dec_factor);


/* wrapper functions */
void interpolate_model(unique_model *in, unique_model *out, interp2d_t *interp_struct);
void gauss_conv_model(unique_model *in, unique_model *out, double stdev);
void decimate_model(unique_model *inv, unique_model *reg, int dec_factor);


#endif