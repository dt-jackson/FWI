#ifndef IBOOL_MODEL_H
#define IBOOL_MODEL_H


typedef struct ibool_model
{
    /* flags */
    int vpvs;
    int kernel;
    int natlog;

    int NPROC;      /* number of processors             */
    int *points;    /* number of points per processor   */
    int *ibool;     /* ibool array                      */
    float *x;       /* x coordinates                    */
    float *z;       /* z coordinates                    */
    float *rho;     /* density                          */
    float *epar1;
    float *epar2; 
    
    float *vp;      /* p-wave velocities                */
    float *vs;      /* s-wave velocities                */

    float *mu;
    float *kappa;

} ibool_model;


void allocate_ibool_model(ibool_model *imod,int NPROC,int *points_per_proc);
void free_ibool_model(ibool_model *imod);

void copy_ibool_model(ibool_model *dest,ibool_model *src);

void read_ibool_model_vpvs(ibool_model *imod,char *path, int NPROC);
void write_ibool_model_vpvs(ibool_model *imod, char *path);

#endif