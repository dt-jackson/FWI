#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


#include <openssl/md5.h> /* for md5 hashing */
/* need to link with libcrypto for md5 */
/* with gcc on linux, add -lcrypto */

#include "fwi.h"
#include "call_specfem.h"
#include "seis_data.h"
#include "unique_model.h"
#include "waveform.h"
#include "traveltime.h"
#include "read_su.h"
#include "write_su.h"
#include "read_kernel.h"
#include "process_seis.h"

#include "mask_kernels.h"
#include "median_filt.h"

#include "fwi_misfit.h"

#define MISFIT_EVAL_MAX 50
#define GRAD_EVAL_MAX 20
#define HASH_SIZE 16   /* md5 hash is 16 bytes */ 

#define MAX_FILE_LEN 256 /* max number of character for file name */

void hash_model(const double *x, const int N, unsigned char *mhash)
{

    /* get size of each element in bytes */
    int ele_bytes = sizeof(*x);

    /* calculate total number of bytes in model */
    long model_bytes = N * ele_bytes;

    /* get pointer to model data */
    unsigned char *modptr = (unsigned char*)x;
    
    /* calculate md5 hash of model */ 
    MD5(modptr, model_bytes, mhash);

    return;
}




double fwi_misfit(const double *x, const int N, void *model)
{
    /* cast model pointer to fwi pointer */
    fwi *ff = (fwi*)model;




    /****************************************************/
    /* hashing stuff */

    static int fevals = 0;
    static unsigned char *fhash = NULL;
    static double *fval = NULL;
    static unsigned char *mhash = NULL;
    //static seis_data *old_resid = NULL;
    //static seis_data old_resid[MISFIT_EVAL_MAX];
    static double *gradient = NULL;

    if (!fhash)
    {
        fhash = fwi_misfit_allocate(sizeof(*fhash) * HASH_SIZE * MISFIT_EVAL_MAX);
    }

    if (!fval)
    {
        fval = fwi_misfit_allocate(sizeof(*fval) * MISFIT_EVAL_MAX);
    }

    //if (!old_resid)
    //{
    //    old_resid = fwi_misfit_allocate(old_resid, MISFIT_EVAL_MAX);
    //}    

    if (!mhash)
    {
        mhash = fwi_misfit_allocate(sizeof(*mhash) * HASH_SIZE);
    }

    if (!gradient)
    {
        gradient = fwi_misfit_allocate(sizeof(*gradient) * N);
    }

    /* hash input model */
    hash_model(x, N, mhash);

    /* write hex hash to string */
    /* each byte will be 2 hex characters, add 1 for trailing null character */
    char mhash_name[2 * HASH_SIZE + 1];
    int kk;
    for (kk=0;kk<HASH_SIZE;kk++)
    {
        sprintf(mhash_name + 2*kk, "%x", mhash[kk]);
    }   

    /* check if model exists already */
    int i;
    int llmin = (MISFIT_EVAL_MAX < fevals) ? MISFIT_EVAL_MAX : fevals; 
    for (i=0;i<llmin;i++)
    {
        if(memcmp(fhash+i*HASH_SIZE, mhash, HASH_SIZE) == 0)
        {
            /* misfit hash been calculated already */

            /* 
            to calculate kernels, we need the correct model, the correct adjoint sources,
            and the following files (in OUTPUT_FILES for each run):
            - lastframe_elastic******.bin
            - pml_interface_elastic******.bin
            - proc******_save_frame_at******.bin
            */

            return fval[i];
        }
    }

    /* done with hashing stuff for now */
    /**************************************************/

    /* evaluate the misfit */
    double f = 0.0;
    
    /* write updated model to bins */
    write_new_model(ff, x);

    /* call specfem to run forward simulation */
    specfem2d_forward(ff->p->sf_bin_forward);

    /* read synthetic seismograms */
    read_synth(ff->syn);

    /* (optional) process synthetic seismograms */
    sr_dist_mute(ff->syn, ff->p->min_sr_dist, ff->p->max_sr_dist);
    process_synth(ff->syn);

    if (ff->p->misfit_type == 2)
    {
        /* traveltime misfit */
        
        /* allocate space for shifts and calculate misfit */
        double *shifts = fwi_misfit_allocate(sizeof(*shifts) * ff->p->num_rec * ff->p->num_src);
        f = traveltime_misfit(ff->obs, ff->syn, shifts);

        /* calcualte adjoint sources */
        traveltime_adjoint(ff->syn, ff->res, shifts);

        /* free shifts */ 
        free(shifts);
    }
    else if (ff->p->misfit_type == 3)
    {
        /* envelope difference misfit */
        f = env_diff_misfit(ff->syn, ff->obs, ff->res);
    }
    else 
    {
        /* waveform misfit -- default case */
    
        /* calculate residual */
        calc_residual(ff->obs, ff->syn, ff->res);
        /* calculate waveform misfit */
        f = waveform_misfit(ff->res, ff->p->data_stdev);
    }


    //print_seis(stderr, ff->res);
    
    /* write adjoint sources */
    write_adjoint_sources(ff->res, ff->p->data_stdev);


    /* call gradient calculation */
    /* we don't need the result here, just need to calculate
    it while the model, adjoint sources, etc. are current */
    fwi_gradient(x, gradient, N, ff);

    /********************************/
    /* hashing stuff */

    /* copy the model hash to the hash array */
    int loc = fevals % MISFIT_EVAL_MAX;
    memcpy(fhash + HASH_SIZE*loc, mhash, HASH_SIZE);

    /* save the adjoint sources */
    //copy_seis_data(old_resid + loc, ff->res);

    /* save the function value to fval */
    fval[loc] = f;

    /* increment fevals and return f */
    fevals++;
    return f;
}





void fwi_gradient(const double *x, double *gradient, const int N, void *model)
{
    fwi *ff = (fwi*)model;
    int ndim = N; 


    /**************************************************/
    /* hashing stuff */

    static int gevals = 0;
    static unsigned char *ghash = NULL;
    static double *gval = NULL;
    static unsigned char *mhash = NULL;


    if (!ghash)
    {
        ghash = fwi_misfit_allocate(sizeof(*ghash) * HASH_SIZE * GRAD_EVAL_MAX);
    }

    if (!gval)
    {
        gval = fwi_misfit_allocate(sizeof(*gval) * ndim * GRAD_EVAL_MAX);
    }

    if (!mhash)
    {
        mhash = fwi_misfit_allocate(sizeof(*mhash) * HASH_SIZE);
    }

    /* hash input model */
    hash_model(x, N, mhash);

    /* check if model exists already */
    int i, j;
    int llmin = GRAD_EVAL_MAX < gevals ? GRAD_EVAL_MAX : gevals; 
    for (i=0;i<llmin;i++)
    {
        if(memcmp(ghash+i*HASH_SIZE, mhash, HASH_SIZE) == 0)
        {
            for (j=0;j<ndim;j++)
            {
                gradient[j] = gval[j + i*ndim];
            }
            return;
        }
    }

    /* done with hashing stuff for now */
    /**************************************************************/

    /* evaluate the gradient */
    /* model and adjoint source should already be written */
    
    /* call specfem to run adjoint simulation */
    specfem2d_adjoint(ff->p->sf_bin_adjoint);

    /* loop through events to read, filter, decimate, and sum kernels */
    char run_no[12];
    char fname[MAX_FILE_LEN];
    int ii;
    for (ii=0;ii<ff->p->num_src;ii++)
    {
        /* write even"t kernel path */
        fname[0] = '\0';
        sprintf(run_no,"%04d",ii+1);
        strcat(fname,"run");
        strcat(fname,run_no);
        strcat(fname,"/OUTPUT_FILES/");

        /* read event kernels */
        read_event_kernel_rmk(ff->imod, fname, ff->p->nproc);        

        /* convert to unique model */
        ibool2unique(ff->imod,ff->umod,ff->i2u);

        double src_loc_x = ff->syn->src_loc[ii];
        double src_loc_z = 0.0;
        double dist_co = ff->p->min_sr_dist;
        double exponent = 2.0;
        mask_src_locs(ff->umod, src_loc_x, src_loc_z, dist_co, exponent);

        /* interpolate to regular grid */
        interpolate_model(ff->umod, ff->umod_reg, ff->sim2reg);

        /* apply median filter */
        int nx = ff->umod_reg->nx;
        int nz = ff->umod_reg->nz;
        int fsize = ff->p->med_filt_size;
        median_filter2d(ff->umod_reg->epar1, nx, nz, fsize, fsize);
        median_filter2d(ff->umod_reg->epar2, nx, nz, fsize, fsize);
        median_filter2d(ff->umod_reg->rho,   nx, nz, fsize, fsize);
        
        /* filter kernels */
        gauss_conv_model(ff->umod_reg, ff->umod_filt, ff->p->sp_lp_cutoff);

        /* decimate kernels */
        int dec_factor = (int)(ff->p->dx_inv / ff->p->dx_reg);
        decimate_model(ff->kernel_inv, ff->umod_filt, dec_factor);

    
        /* mask kernels */
        if (ff->p->mask_dens)
            mask_density(ff->kernel_inv);
        if (ff->p->mask_kernels)
            mask_kernel(ff->kernel_inv, ff->mask);
        



        /* sum event kernels - put result in gradient vector */
        if (ii==0)
        {
            for (j=0;j<ndim;j++)
            {
                gradient[j] = ff->kernel_inv->data[j];
            }
        }
        else
        {
            for (j=0;j<ndim;j++)
            {
                gradient[j] += ff->kernel_inv->data[j];
            }
        }

    }
     

    /*******************************************************************/
    /* mode hashing stuff here */

    /* copy the model hash to the hash array */
    int loc = gevals % GRAD_EVAL_MAX;
    memcpy(ghash + HASH_SIZE*loc, mhash, HASH_SIZE);

    /* save the gradient value to gval */
    for (j=0;j<ndim;j++)
    {
        gval[j + loc*ndim] = gradient[j];
    }

    /* increment gevals and return */
    gevals++;

    return;
}

/*****************************************************************************/



int get_ndim(void *model)
{
    /* multiply by 3 because there are 3 parameters */
    /* density, vp, vs or density, mu, kappa        */

    return 3 * ((fwi*)model)->umod_inv->num_ele;
}

void model2vector(void *model, double *x)
{
    fwi *ff = (fwi*)model;
    
    /* multiply by 3 because there are 3 parameters */
    /* density, vp, vs or density, mu, kappa        */
    int N = 3 * (ff->umod_inv->num_ele);

    int i;
    for (i=0;i<N;i++)
    {
        x[i] = ff->umod_inv->data[i];
    }

}


void vector2model(void *model, const double *x)
{
    fwi *ff = (fwi*)model;
    
    /* multiply by 3 because there are 3 parameters */
    /* density, vp, vs or density, mu, kappa        */
    int N = 3* (ff->umod_inv->num_ele);

    int i;
    for (i=0;i<N;i++)
    {
        ff->umod_inv->data[i] = x[i];
    }

}


void write_adjoint_sources(seis_data *d, double stdev)
{
    /* scale adjoint sources by inverse variance */
    double inv_variance = 1.0 / (stdev * stdev);
    seis_data tmp;
    tmp.src_loc = NULL;
    tmp.rec_loc = NULL;
    tmp.traces = NULL;
    copy_seis_data(&tmp, d);

    int total_samp = tmp.num_samples * tmp.num_rec * tmp.num_src;
    int j;
    for (j=0;j<total_samp;j++)
    {
        tmp.traces[j] *= inv_variance;
    }

    int i;
    char src_no[12];
    char fname[MAX_FILE_LEN];
    for (i=0;i<d->num_src;i++)
    {
        fname[0] = '\0';
        sprintf(src_no,"%04d",i+1);
        strcat(fname,"run");
        strcat(fname,src_no);
        strcat(fname,"/SEM/Uz_file_single.su.adj");
        write_su_data_only(fname, &tmp, i);
    }

    free_seis_data(&tmp);
}

void write_new_model(fwi *ff, const double *x)
{

    vector2model(ff,x);

    /* convert model to sim grid */  
    interpolate_model(ff->umod_inv, ff->umod, ff->inv2sim);

    /* convert to vpvs if necessary */
    exp_mod(ff->umod); 
    rmk2vpvs(ff->umod);

    enforce_par_bounds(ff->umod, ff->p);

    /* convert model to ibool */
    unique2ibool(ff->imod, ff->umod, ff->i2u);

    /* write ibool model to bins */
    write_ibool_model_vpvs(ff->imod, ff->p->model_path);

}


void read_synth(seis_data *s)
{
    int i;
    char src_no[12];
    char fname[MAX_FILE_LEN];
    for (i=0;i<s->num_src;i++)
    {
        fname[0] = '\0';
        sprintf(src_no,"%04d",i+1);
        strcat(fname,"run");
        strcat(fname,src_no);
        strcat(fname,"/OUTPUT_FILES/Uz_file_single.su");
        read_su_data_only(fname, s, i);
    }

}


void process_synth(seis_data *s)
{
    /* placeholder function for now */
    return;
}


void* fwi_misfit_allocate(long nbytes)
{
    void *ptr = malloc(nbytes);
    if (!ptr)
    {
        fprintf(stderr,"Error: fwi_misfit: memory allocation failure\n");
        exit(1);
    }
    return ptr;
}


void enforce_par_bounds(unique_model *mod, parameters *p)
{
    /* ad hoc enforcement of parameter space bounds */
    double vp_vs_ratio_max = sqrt(2*(1-p->max_pr) / (1-2*p->max_pr));
    double vp_vs_ratio_min = sqrt(2*(1-p->min_pr) / (1-2*p->min_pr));
    double vp_vs_ratio;
 
    int k;
    for (k=0;k<mod->num_ele;k++)
    {

        /* enforce vp bounds */
        if (mod->vp[k] > p->max_vp)
            mod->vp[k] = p->max_vp;

        if (mod->vp[k] < p->min_vp)
            mod->vp[k] = p->min_vp;

        /* enforce vs bounds */
        if (mod->vs[k] > p->max_vs)
            mod->vs[k] = p->max_vs;

        if (mod->vs[k] < p->min_vs)
            mod->vs[k] = p->min_vs;
        

        /* enforce poisson ratio bounds */
        vp_vs_ratio = mod->vp[k] / mod->vs[k];

        if (vp_vs_ratio > vp_vs_ratio_max)
        {
            /* set poisson equal to max_pr, fix vs, and solve for vp */
            mod->vp[k] = mod->vs[k] * vp_vs_ratio_max;
        } 
        if (vp_vs_ratio < vp_vs_ratio_min)
        {
            /* set poisson equal to min_pr, fix vs, and solve for vp */
            mod->vp[k] = mod->vs[k] * vp_vs_ratio_min;
        }


        /* enforce rho bounds */            
        if (mod->rho[k] > p->max_rho)
            mod->rho[k] = p->max_rho;

        if (mod->rho[k] < p->min_rho)
            mod->rho[k] = p->min_rho;
    }



}



void mask_src_locs(unique_model *kernel, double src_loc_x, double src_loc_z, double dist_co, double exponent)
{

    int  i;
    int num_ele = kernel->num_ele;
    
    
    double dist_co2 = dist_co * dist_co;


    double dx, dz, dist, dist2, mask;
    for (i=0;i<num_ele;i++)
    {

        dx = (kernel->x[i] - src_loc_x);
        dz = (kernel->z[i] - src_loc_z);
        dist2 = dx*dx + dz*dz;


        if (dist2 < dist_co2)
        {
            dist = sqrt(dist2);
            mask = pow((dist / dist_co), exponent);

            kernel->epar1[i] *= mask;
            kernel->epar2[i] *= mask;
            kernel->dens[i] *= mask;
        }
    }





    return;
}