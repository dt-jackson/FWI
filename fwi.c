#include <stdlib.h>
#include <stdio.h>

#include "fwi.h"
#include "seis_data.h"
#include "unique_model.h"
#include "bilinear_interp.h"
#include "gauss_conv.h"

#include "mask_kernels.h"

#include "write_su.h"
#include "read_su.h"

#include "read_src_rec_loc.h"

#include "process_seis.h"

#define MAX_FILE_LEN 256

/*****************************************************************************/

void initialize_fwi(fwi *ff, parameters *p)
{
    ff->p = p;

    /* allocate observed and synthetic data */
    ff->obs = &ff->obs_;
    ff->syn = &ff->syn_;
    ff->res = &ff->res_;
    allocate_seis_data(ff->obs, p->num_src, p->num_rec, p->nstep);
    allocate_seis_data(ff->syn, p->num_src, p->num_rec, p->nstep);
    allocate_seis_data(ff->res, p->num_src, p->num_rec, p->nstep);
    ff->obs->dt = p->dt;
    ff->syn->dt = p->dt;
    ff->res->dt = p->dt;

    read_src_loc(p, ff->obs);
    read_rec_loc(p, ff->obs);
    int kk;
    for (kk=0;kk<ff->obs->num_src;kk++)
    {
        ff->syn->src_loc[kk] = ff->obs->src_loc[kk];
        ff->res->src_loc[kk] = ff->obs->src_loc[kk];
    }
    for (kk=0;kk<ff->obs->num_rec;kk++)
    {
        ff->syn->rec_loc[kk] = ff->obs->rec_loc[kk];
        ff->res->rec_loc[kk] = ff->obs->rec_loc[kk];
    }



    //ff->res->src_loc[0] = 5;
    //ff->res->src_loc[1] = 20;


    //for (int i=0;i<27;i++)
    //{
    //    ff->res->rec_loc[i] = (double)i;
    //}
    

    /* read observed data */
    char src_no[12];
    char fname[MAX_FILE_LEN];
    int ii;
    for (ii=0;ii<ff->obs->num_src;ii++)
    {
        fname[0] = '\0';
        sprintf(src_no,"%04d",ii+1);
        strcat(fname,"run");
        strcat(fname,src_no);
        strcat(fname,"/Uz_file_single.su");
        read_su_data_only(fname, ff->obs, ii);
    }

    /* mute near and far traces */
    sr_dist_mute(ff->obs, p->min_sr_dist, p->max_sr_dist);





    /* write dummy adjoint source files */
    for (ii=0;ii<ff->syn->num_src;ii++)
    {
        fname[0] = '\0';
        sprintf(src_no,"%04d",ii+1);
        strcat(fname,"run");
        strcat(fname,src_no);
        strcat(fname,"/SEM/Ux_file_single.su.adj");
        write_su_data_only(fname,ff->res,-1);

        fname[0] = '\0';
        strcat(fname,"run");
        strcat(fname,src_no);
        strcat(fname,"/SEM/Uy_file_single.su.adj");
        write_su_data_only(fname,ff->res,-1);
    }


    /* read initial model - imod will be allocated here*/
    ff->imod = &ff->imod_;
    read_ibool_model_vpvs(ff->imod, p->model_path, p->nproc);

    /* construct conversion between ibool and unique models */
    ff->i2u = &ff->i2u_;
    construct_ibool_unique(ff->imod,ff->i2u);

    /* allocate unique model */
    /* conversion is performed here so coordinates can be   */
    /* used for interpolation                               */
    ff->umod = &ff->umod_;
    allocate_unique_model(ff->umod, ff->i2u->num_unique); 
    ibool2unique(ff->imod, ff->umod, ff->i2u);

    /* get grid sizes */
    int nx_reg, nz_reg, nx_inv, nz_inv;
    get_grid_size(p, &nx_reg, &nz_reg, &nx_inv, &nz_inv);

    /* make regular grid model */
    ff->umod_reg = &ff->umod_reg_;
    allocate_unique_model(ff->umod_reg, nx_reg*nz_reg);
    make_reg_grid(ff->umod_reg,p->dx_reg, p->xmin, p->xmax, p->zmin, p->zmax); 
    ff->umod_reg->nx = nx_reg;
    ff->umod_reg->nz = nz_reg;
   

    /* allocate filtered grid */
    ff->umod_filt = &ff->umod_filt_;
    allocate_unique_model(ff->umod_filt, ff->umod_reg->num_ele);
    ff->umod_filt->nx = nx_reg;
    ff->umod_filt->nz = nz_reg;

    /* make decimated grid -allocatd here */
    ff->umod_inv = &ff->umod_inv_;
    int dec_factor = (int)(p->dx_inv / p->dx_reg);
    allocate_unique_model(ff->umod_inv, nx_inv * nz_inv);
    decimate_grid(ff->umod_inv, ff->umod_reg,dec_factor);
    ff->umod_inv->nx = nx_inv;
    ff->umod_inv->nz = nz_inv;


    /* allocate decimated kernels */
    ff->kernel_inv = &ff->kernel_inv_;
    allocate_unique_model(ff->kernel_inv, ff->umod_inv->num_ele);
    ff->kernel_inv->nx = nx_inv;
    ff->kernel_inv->nz = nz_inv;
    

    /* make interpolation from unique model to regular grid */
    ff->sim2reg = malloc(sizeof(*ff->sim2reg) * (ff->umod_reg->num_ele));
    if (!ff->sim2reg)
    {
        fprintf(stderr,"Memory allocation failure: sim2reg\n");
        exit(-1);
    }
    bilinear_setup(ff->umod->x, ff->umod->z, ff->umod_reg->x, ff->umod_reg->z,
                   ff->umod->num_ele, ff->umod_reg->num_ele, ff->sim2reg);

    /* make interpolation from inversion grid to simulation grid */
    ff->inv2sim = malloc(sizeof(*ff->inv2sim)*ff->umod->num_ele);
    if (!ff->inv2sim)
    {
        fprintf(stderr,"Memory allocation failure: inv2sim\n");
        exit(-1);
    }
    bilinear_setup(ff->umod_inv->x, ff->umod_inv->z, ff->umod->x, ff->umod->z,
                   ff->umod_inv->num_ele, ff->umod->num_ele, ff->inv2sim);


    /* set up mask */
    ff->mask = &ff->mask_;
    allocate_unique_model(ff->mask, ff->umod_inv->num_ele);
    if (p->mask_kernels)
        prepare_mask(p->mask_xmin, p->mask_xmax, p->mask_zmin, p->mask_zmax, ff->umod_inv, ff->mask);

}

/*****************************************************************************/
void free_fwi(fwi *ff)
{
    if (!ff)
        return;

    if (ff->obs)
        free_seis_data(ff->obs);
    if (ff->syn) 
        free_seis_data(ff->syn);
    if (ff->res) 
        free_seis_data(ff->res);

    if (ff->imod) 
        free_ibool_model(ff->imod);
    if (ff->i2u)
        free_conversion_struct(ff->i2u);
    if (ff->umod)
        free_unique_model(ff->umod);

    if (ff->umod_reg)
        free_unique_model(ff->umod_reg);
    if (ff->umod_filt)
        free_unique_model(ff->umod_filt);

    if (ff->umod_inv) 
        free_unique_model(ff->umod_inv);
    if (ff->kernel_inv)
        free_unique_model(ff->kernel_inv);

    if (ff->sim2reg)
        free(ff->sim2reg);
    if (ff->inv2sim)
        free(ff->inv2sim);

    if (ff->mask)
        free_unique_model(ff->mask);
}

/*****************************************************************************/

int get_dec_factor(parameters *p)
{
    return (int)(p->dx_inv / p->dx_reg);
}


void get_grid_size(parameters *p,int *nx_reg, int *nz_reg, int *nx_inv, int *nz_inv)
{
    int dec_factor = (int)(p->dx_inv / p->dx_reg);

    /* regular grid size */
    *nx_reg = ((p->xmax - p->xmin) / p->dx_reg) + 1;
    *nz_reg = ((p->zmax - p->zmin) / p->dx_reg) + 1;

    /* inversion grid size */
    *nx_inv = (*nx_reg / dec_factor) + 1;
    *nz_inv = (*nz_reg / dec_factor) + 1;
}

void make_reg_grid(unique_model *reg, double dx, double xmin, double xmax, double zmin, double zmax)
{

    /* get number of points in each direction */
    int nx = ((xmax - xmin) / dx) + 1;
    int nz = ((zmax - zmin) / dx) + 1;

    int total_points = nx * nz;


    /* write the points to the model */ 
    /* we're only going through this loop once, so the branches
    shouldn't be a big deal */ 
    int i, j, index;
    double x, z;
    for (j=0;j<nz;j++)
    {
        z = zmin + (dx * j);

        /* assign max and min values to make sure
        there is no error on the boundary due to 
        floating point arithemetic */
        if (j==0)
            z = zmin;
        if (j==nz-1)
            z = zmax;

        for (i=0;i<nx;i++)
        {
            x = xmin + (dx * i);

            if (i==0)
                x = xmin;
            if (i==nx-1)
                x = xmax;

            index = i + (j * nx);
            reg->x[index] = x;
            reg->z[index] = z;
        }
    }
    reg->nx = nx;
    reg->nz = nz;
}


/*****************************************************************************/
void decimate_grid(unique_model *inv, unique_model *reg, int dec_factor)
{
    int nx = (reg->nx / dec_factor) + 1;
    int nz = (reg->nz / dec_factor) + 1;

    int i, j, index_inv, index_reg;
    for (j=0;j<nz;j++)
    {
        for (i=0;i<nx;i++)
        {
            index_inv = i + (j * nx);
            index_reg = (dec_factor * i) + (dec_factor * j * reg->nx);
            inv->x[index_inv] = reg->x[index_reg];
            inv->z[index_inv] = reg->z[index_reg];
        }
    }
    inv->nx = nx;
    inv->nz = nz;

}






/*****************************************************************************/
/* wrapper functions */

void interpolate_model(unique_model *in, unique_model *out, interp2d_t *interp_struct)
{
    bilinear_interp(in->epar1, out->epar1, out->num_ele, interp_struct);
    bilinear_interp(in->epar2, out->epar2, out->num_ele, interp_struct);
    bilinear_interp(in->dens,  out->dens,  out->num_ele, interp_struct);

    /* copy flags */ 
    out->vpvs   = in->vpvs;
    out->natlog = in->natlog;
    out->kernel = in->kernel;
}

void gauss_conv_model(unique_model *in, unique_model *out, double stdev)
{
    if (in->nx <= 0 || in->nz <= 0)
    {
        fprintf(stderr,"Grid size not assigned\n");
        exit(1);
    }

    gauss_conv2d(in->epar1, out->epar1, in->nx, in->nz, stdev, stdev);
    gauss_conv2d(in->epar2, out->epar2, in->nx, in->nz, stdev, stdev);
    gauss_conv2d(in->dens,  out->dens,  in->nx, in->nz, stdev, stdev);

    /* copy flags */ 
    out->vpvs   = in->vpvs;
    out->natlog = in->natlog;
    out->kernel  = in->kernel;

    /* copy sizes */
    out->num_ele = in->num_ele;
    out->nx = in->nx;
    out->nz = in->nz;
}


void decimate_model(unique_model *inv, unique_model *reg, int dec_factor)
{
    int nx = inv->nx;
    int nz = inv->nz;

    int i, j, index_inv, index_reg;
    for (j=0;j<nz;j++)
    {
        for (i=0;i<nx;i++)
        {
            index_inv = i + (j * nx);
            index_reg = (dec_factor * i) + (dec_factor * j * reg->nx);
            inv->epar1[index_inv] = reg->epar1[index_reg];
            inv->epar2[index_inv] = reg->epar2[index_reg];
            inv->dens[index_inv]  = reg->dens[index_reg];
        }
    }

    /* copy flags */
    inv->vpvs   = reg->vpvs;
    inv->natlog = reg->natlog;
    inv->kernel = reg->kernel;

}

