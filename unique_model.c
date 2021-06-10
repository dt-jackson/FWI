#include <stdlib.h>
#include <stdio.h>
#include <math.h>           /*for fabs, sqrt, log, exp*/

#include "ibool_model.h"

#include "unique_model.h"


#define TOL 1e-6 /* tolerance beteen points when checking uniqueness */

/*****************************************************************************/
void vpvs2rmk(unique_model *umod)
{
    if (umod->kernel)
    {
        fprintf(stderr,"vp/vs to rmk conversion should not be performed on kernels\n");
        exit(1);
    }
    if (!umod->vpvs)
    {
        //fprintf(stderr,"Model not vp/vs\n");
        //exit(-50);
        return;
    }

    if (umod->natlog)
    {
        exp_mod(umod);
    }


    long i;
    long N = umod->num_ele;

    double rho;
    double M, mu, kappa;
    double vp, vs;
    double four_thirds = 4.0 / 3.0; 
    for (i=0;i<N;i++)
    {
        rho = umod->dens[i];
        vp = umod->epar1[i]; 
        vs = umod->epar2[i];
        mu = vs * vs * rho;
        M = vp * vp * rho;
        kappa = M - four_thirds * mu;

        umod->mu[i] = mu;
        umod->kappa[i] = kappa;
    }
    umod->vpvs = 0;
}

void rmk2vpvs(unique_model *umod)
{
    if (umod->kernel)
    {
        fprintf(stderr,"rmk to vp/vs conversion should not be performed on kernels\n");
        exit(1);
    }

    if (umod->vpvs)
    {
        //fprintf(stderr,"Model not rmk\n");
        //exit(-50);
        return;
    }

    if (umod->natlog)
    {
        exp_mod(umod);
    }


    long i;
    long N = umod->num_ele;

    double rho; 
    double M, mu, kappa;
    double vp, vs;
    double four_thirds = 4.0 / 3.0; 

    for (i=0;i<N;i++)
    {
        rho = umod->rho[i];
        mu = umod->mu[i];
        kappa = umod->kappa[i];
        M = kappa + four_thirds * mu;

        umod->vp[i] = sqrt(M / rho);
        umod->vs[i] = sqrt(mu / rho);
    }
    umod->vpvs = 1;
}

void ln_mod(unique_model *umod)
{
    if (umod->kernel)
    {
        fprintf(stderr,"natural log conversion should not be performed on kernels\n");
        exit(1);
    }
    if (umod->natlog)
    {
        //fprintf(stderr,"Model already natural log\n");
        //exit(-50);
        return;
    }

    long i;
    long N = umod->num_ele;

    for (i=0;i<N;i++)
    {
        umod->dens[i] = log(umod->dens[i]);
        umod->epar1[i] = log(umod->epar1[i]);
        umod->epar2[i] = log(umod->epar2[i]);
    }
    umod->natlog = 1;
}

void exp_mod(unique_model *umod)
{
    if (umod->kernel)
    {
        fprintf(stderr,"exp conversion should not be performed on kernels\n");
        exit(1);
    }
    if (!umod->natlog)
    {
        //fprintf(stderr,"Model not natural log\n");
        //exit(-50);
        return;
    }
    
    long i;
    long N = umod->num_ele;

    for (i=0;i<N;i++)
    {
        umod->dens[i] = exp(umod->dens[i]);
        umod->epar1[i] = exp(umod->epar1[i]);
        umod->epar2[i] = exp(umod->epar2[i]);
    }
    umod->natlog = 0;
}



/*****************************************************************************/
/* Functions for performing conversion                                       */

void ibool2unique(ibool_model *imod, unique_model *umod,
                  convert_ibool_unique *conversion)
{
    /* Note that umod must be allocated already */
    int i;
    for (i=0;i<conversion->num_unique;i++)
    {
        umod->x[i]   = (double)imod->x[conversion->ibool2unique[i]];
        umod->z[i]   = (double)imod->z[conversion->ibool2unique[i]];
        umod->epar1[i]  = (double)imod->epar1[conversion->ibool2unique[i]];
        umod->epar2[i]  = (double)imod->epar2[conversion->ibool2unique[i]];
        umod->rho[i] = (double)imod->rho[conversion->ibool2unique[i]];
    }
    umod->vpvs = imod->vpvs;
    umod->kernel = imod->kernel;
    umod->natlog = imod->natlog;
}

void unique2ibool(ibool_model *imod, unique_model *umod,
                  convert_ibool_unique *conversion)
{
    /* Note that imod must be allocated already */
    int i;
    for (i=0;i<conversion->num_ibool;i++)
    {
        imod->x[i]   = (float)umod->x[conversion->unique2ibool[i]];
        imod->z[i]   = (float)umod->z[conversion->unique2ibool[i]];
        imod->epar1[i]  = (float)umod->epar1[conversion->unique2ibool[i]];
        imod->epar2[i]  = (float)umod->epar2[conversion->unique2ibool[i]];
        imod->rho[i] = (float)umod->rho[conversion->unique2ibool[i]];
    }
    imod->vpvs = umod->vpvs;
    imod->kernel = umod->kernel;
    imod->natlog = umod->natlog;
}




/*****************************************************************************/
/* Functions for cconstructing conversion                                    */

void construct_ibool_unique(ibool_model *imod,
                           convert_ibool_unique *conversion)
{
    int i;

    /* get total number of elements in ibool model */    
    int no_ele = 0;
    for (i=0;i<imod->NPROC;i++)
    {
        no_ele += imod->points[i];
    }

    /* copy coordinaates into sortable struct with index */
    coord *ind = calloc(sizeof(*ind), no_ele);
    if (!ind)
    {
        fprintf(stderr,"Memory allocation failure: construct_ibool_unique: ind\n");
        exit(1);
    }

    for (i=0;i<no_ele;i++)
    {
        ind[i].index = i;
        ind[i].x = imod->x[i];
        ind[i].z = imod->z[i];
    } 

    /* sort coordinates */
    qsort(ind,no_ele,sizeof(*ind),sort_by_coord);

    /* unique2ibool converts from an array of unique values back */
    /* to the full ibool values                                  */
    int *u2i = calloc(no_ele, sizeof(int));
    if (!u2i)
    {
        fprintf(stderr,"Memory allocation failure: construct_ibool_unique: u2i\n");
        exit(1);
    }
    u2i[ind[0].index] = 0;

    /* get index of unique elements */
    int count = 1;
    for (i=1;i<no_ele;i++)
    {
        if (compare_coord(&ind[i-1],&ind[i]) == 0)
        {
            ind[count++] = ind[i];
        }
        u2i[ind[i].index] = count - 1;
    }

    /* allocate output struct and write values */
    allocate_conversion_struct(conversion, count, no_ele);

    for (i=0;i<count;i++)
    {
        conversion->ibool2unique[i] = ind[i].index;
    }

    for (i=0;i<no_ele;i++)
    {
        conversion->unique2ibool[i] = u2i[i];
    }


    free(u2i);
    free(ind);

}


int sort_by_coord(const void *p1, const void *p2)
{
    /* returns -1 if p1 goes before p2 */
    /* returns  1 if p1 goes after p2  */

    double x1 = ((coord*)p1)->x;
    double x2 = ((coord*)p2)->x;

    double z1 = ((coord*)p1)->z;
    double z2 = ((coord*)p2)->z;

    /* If z1 = z2, sort by x. Otherwise, sort by z */
    if (fabs(z1 - z2) < TOL)
    {
        /* z values are equal */
        if (x1 < x2)
            return -1;
        else
            return 1;
    }
    else
    {
        /* z values are different */
        if (z1 < z2)
            return -1;
        else 
            return 1; 
    }
}

int compare_coord(const void *p1,const void *p2)
{
    /* returns 1 if p1 and p2 are equal */
    /* returns 0 otherwise              */

    if (fabs(((coord*)p1)->x - ((coord*)p2)->x) < TOL 
        && fabs(((coord*)p1)->z -((coord*)p2)->z) < TOL)
    {
        return 1;
    }

    return 0;
}


/*****************************************************************************/
/* Functions for memory management                                           */

void allocate_conversion_struct(convert_ibool_unique *conversion,
                                int num_unique, int num_ibool)
{
    conversion->num_ibool = num_ibool;
    conversion->num_unique = num_unique;

    //conversion->ibool2unique = calloc(num_unique, sizeof(*conversion->ibool2unique));
    //conversion->unique2ibool = calloc(num_ibool, sizeof(*conversion->unique2ibool));
    conversion->ibool2unique = malloc(sizeof(*conversion->ibool2unique) * num_unique);
    conversion->unique2ibool = malloc(sizeof(*conversion->unique2ibool) * num_ibool);
    if (!conversion->ibool2unique || !conversion->unique2ibool)
    {
        fprintf(stderr,"Memory allocation failure: allocate_conversion_struct\n");
        exit(1);
    }
}

void free_conversion_struct(convert_ibool_unique *conversion)
{
    if (!conversion)
        return;
    
    if (conversion->ibool2unique)
    {
        free(conversion->ibool2unique);
        conversion->ibool2unique = NULL;
    }
    if (conversion->unique2ibool)
    {
        free(conversion->unique2ibool);
        conversion->unique2ibool = NULL;
    }
}


void allocate_unique_model(unique_model *umod,int num_ele)
{

    umod->num_ele = num_ele;
    umod->nx   = -1; /* start out with nx = nz = -1 to indicate it has not been assigned */
    umod->nz   = -1;
    umod->x    = malloc(sizeof(*umod->x) * num_ele);
    umod->z    = malloc(sizeof(*umod->z) * num_ele);
    umod->data = malloc(sizeof(*umod->data) * 3 * num_ele);
    /* check if allocation succeeded */
    if (!umod->x || !umod->z || !umod->data)
    {
        fprintf(stderr,"Memory allocation failure: allocate_unique_model\n");
        exit(1);
    }

    umod->epar1 = umod->data;
    umod->epar2 = umod->data + num_ele;
    umod->dens  = umod->data + 2*num_ele;

    umod->vpvs = 0;
    umod->natlog = 0;
    umod->kernel = 0;
    
    umod->rho = umod->dens;

    umod->mu = umod->epar1;
    umod->kappa = umod->epar2;

    umod->vp = umod->epar1;
    umod->vs = umod->epar2;
}

void free_unique_model(unique_model *umod)
{
    if (!umod)
        return;
    
    if (umod->x)
    {
        free(umod->x);
        umod->x = NULL;
    }

    if (umod->z)
    {
        free(umod->z);
        umod->z = NULL;
    }

    if (umod->data)
    {
        free(umod->data);
        umod->data = NULL;
    }
}

/*****************************************************************************/

void copy_unique_model(unique_model *dest, unique_model *src)
{
    int N = src->num_ele;
    
    /* free and reallocate dest model */
    free_unique_model(dest);
    allocate_unique_model(dest, N);

    /* copy sizes */
    dest->num_ele = src->num_ele;
    dest->nx = src->nx;
    dest->nz = src->nz;


    /* copy flags */
    dest->vpvs = src->vpvs;
    dest->natlog = src->natlog;
    dest->kernel = src->kernel;

    /* deep copy of arrays */
    int i;
    for (i=0;i<N;i++)
    {
        dest->x[i] = src->x[i]; 
        dest->z[i] = src->z[i]; 
        dest->rho[i] = src->rho[i]; 
        dest->epar1[i] = src->epar1[i]; 
        dest->epar2[i] = src->epar2[i]; 
    }
    



}







/*****************************************************************************/


void model_min_max(unique_model *umod)
{
    double dens_min, dens_max;
    double epar1_min, epar1_max;
    double epar2_min, epar2_max;


    dens_min = umod->dens[0];
    dens_max = umod->dens[0];

    epar1_min = umod->epar1[0];
    epar1_max = umod->epar1[0];

    epar2_min = umod->epar2[0];
    epar2_max = umod->epar2[0];

    int N = umod->num_ele;
    int i;
    for (i=0;i<N;i++)
    {
        dens_min = dens_min < umod->dens[i] ? dens_min : umod->dens[i];
        dens_max = dens_max > umod->dens[i] ? dens_max : umod->dens[i];

        epar1_min = epar1_min < umod->epar1[i] ? epar1_min : umod->epar1[i];
        epar1_max = epar1_max > umod->epar1[i] ? epar1_max : umod->epar1[i];

        epar2_min = epar2_min < umod->epar2[i] ? epar2_min : umod->epar2[i];
        epar2_max = epar2_max > umod->epar2[i] ? epar2_max : umod->epar2[i];
    }


    printf("Density min/max : %E %E\n", dens_min, dens_max);
    printf("Par1 min/max    : %E %E\n", epar1_min, epar1_max);
    printf("Par2 min/max    : %E %E\n", epar2_min, epar2_max);

}