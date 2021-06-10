#ifndef UNIQUE_MODEL_H
#define UNIQUE_MODEL_H

#include "ibool_model.h"

typedef struct unique_model
{
    /* flags, 1 if true, 0 if false */
    int vpvs;
    int natlog;
    int kernel;

    /* sizes */
    int num_ele;
    int nx;
    int nz;

    /* data */
    double *x;
    double *z;
    double *data;
    double *epar1;
    double *epar2;
    double *dens;

    /* pointers to data */
    double *rho;    /* equal to dens                 */
    double *mu;     /* either equal to epar1 or NULL */
    double *kappa;  /* either equal to epar2 or NULL */
    double *vp;     /* either equal to epar1 or NULL */
    double *vs;     /* either equal to epar2 or NULL */

} unique_model;




typedef struct convert_ibool_unique
{
    int num_ibool;
    int num_unique;
    int *ibool2unique;
    int *unique2ibool;
} convert_ibool_unique;

typedef struct coord
{
    int index;
    double x;
    double z;
} coord;


void ibool2unique(ibool_model *imod, unique_model *umod,
                  convert_ibool_unique *conversion);

void unique2ibool(ibool_model *imod, unique_model *umod,
                  convert_ibool_unique *conversion);


void construct_ibool_unique(ibool_model *imod,
                           convert_ibool_unique *conversion);
int sort_by_coord(const void *p1, const void *p2);
int compare_coord(const void *p1,const void *p2);


void allocate_conversion_struct(convert_ibool_unique *conversion,
                                int num_unique, int num_ibool);
void free_conversion_struct(convert_ibool_unique *conversion);

void allocate_unique_model(unique_model *umod,int num_ele);
void free_unique_model(unique_model *umod);
void copy_unique_model(unique_model *dest, unique_model *src);

void vpvs2rmk(unique_model *umod);
void rmk2vpvs(unique_model *umod);

void ln_mod(unique_model *umod);
void exp_mod(unique_model *umod);



void model_min_max(unique_model *umod);


#endif