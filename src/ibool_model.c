#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "ibool_model.h"
#include "bin_io.h"


#define MAX_FNAME 512

/*****************************************************************************/
void allocate_ibool_model(ibool_model *imod,int NPROC,int *points_per_proc)
{
    int i;
    
    imod->NPROC = NPROC;
    imod->points = calloc(NPROC, sizeof(int));

    if (imod->points == NULL)
    {
        /*TODO: Handle mem allocation failure*/
    }

    int total_pts = 0;
    for (i=0;i<NPROC;i++)
    {
        imod->points[i] = points_per_proc[i];
        total_pts += points_per_proc[i];
    }

    imod->ibool = calloc(total_pts, sizeof(int));
    imod->x = calloc(total_pts, sizeof(float));
    imod->z = calloc(total_pts, sizeof(float));
    imod->epar1 = calloc(total_pts, sizeof(float));
    imod->epar2 = calloc(total_pts, sizeof(float));
    imod->rho = calloc(total_pts, sizeof(float));
    
    if (imod->x == NULL || imod->z == NULL || imod->epar1 == NULL ||
        imod->epar2 == NULL || imod->rho == NULL)
    {
        /*TODO: Handle mem allocation failure*/
    }

    /* vp/vs by default */
    imod->vpvs = 1;
    imod->vp = imod->epar1;
    imod->vs = imod->epar2;
    imod->mu = imod->epar1;
    imod->kappa = imod->epar2;
}

/*****************************************************************************/

void free_ibool_model(ibool_model *imod)
{
    if (!imod)
        return;

    if (imod->points)
    {
        free(imod->points);
        imod->points = NULL;
    }
    
    if (imod->ibool)
    {
        free(imod->ibool);
        imod->ibool = NULL;
    }

    if (imod->x)
    {
        free(imod->x);
        imod->x = NULL;
    }
    if (imod->z)
    {
        free(imod->z);
        imod->z = NULL;
    }
    
    if (imod->epar1)
    {
        free(imod->epar1);
        imod->epar1 = NULL;
    }
    if (imod->epar2)
    {
        free(imod->epar2);
        imod->epar2 = NULL;
    }
    if (imod->rho)
    {
        free(imod->rho);
        imod->rho = NULL;
    }
}

/*****************************************************************************/

void copy_ibool_model(ibool_model *dest,ibool_model *src)
{
    int i, j;
    int total_pts = 0;

    /* allocate new model. NPROC and points will be copied here */
    free_ibool_model(dest);
    allocate_ibool_model(dest,src->NPROC,src->points);

    /* get total number of elements in each array */
    for (i=0;i<src->NPROC;i++)
    {
        total_pts += src->points[i];
    }

    /* deep copy of arrays */
    for (i=0;i<total_pts;i++)
    {
        dest->ibool[i] = src->ibool[i];
        dest->x[i] = src->x[i];
        dest->z[i] = src->z[i];
        dest->epar1[i] = src->epar1[i];
        dest->epar2[i] = src->epar2[i];
        dest->rho[i] = src->rho[i];
    }

    /* copy flags */
    dest->vpvs = src->vpvs;
    dest->natlog = src->natlog;
    dest->kernel = src->kernel;


}


/*****************************************************************************/

void read_ibool_model_vpvs(ibool_model *imod,char *path, int NPROC)
{
    /* set flags */
    imod->vpvs = 1;
    imod->natlog = 0;
    imod->kernel = 0;

    int i;

    char fname[MAX_FNAME];
    char tmp_string[MAX_FNAME];


    /* initialize arrays of file pointers */
    FILE **x_files     = (FILE**)calloc(NPROC, sizeof(FILE*));
    FILE **z_files     = (FILE**)calloc(NPROC, sizeof(FILE*));
    FILE **ibool_files = (FILE**)calloc(NPROC, sizeof(FILE*));
    FILE **vp_files    = (FILE**)calloc(NPROC, sizeof(FILE*));
    FILE **vs_files    = (FILE**)calloc(NPROC, sizeof(FILE*));
    FILE **rho_files   = (FILE**)calloc(NPROC, sizeof(FILE*));

    /* initialize local array for number of points to be read */
    int *points_per_proc = (int*)calloc(NPROC, sizeof(int));

    /* First loop to open files and get array sizes*/
    for (i=0;i<NPROC;i++)
    {
        /* write processor number to string */
        char proc_no[12];
        sprintf(proc_no,"%06d",i);

        /* write base file name */
        strcpy(fname,path);
        strcat(fname,"/proc");
        strcat(fname,proc_no);

        /* open x file */ 
        strcpy(tmp_string,fname);
        strcat(tmp_string,"_x.bin");
        x_files[i] = open_bin_read(tmp_string);

        /* open z file */
        strcpy(tmp_string,fname);
        strcat(tmp_string,"_z.bin");
        z_files[i] = open_bin_read(tmp_string);

        /* open ibool file */
        strcpy(tmp_string,fname);
        strcat(tmp_string,"_NSPEC_ibool.bin");
        ibool_files[i] = open_bin_read(tmp_string);

        /* open vp file */ 
        strcpy(tmp_string,fname);
        strcat(tmp_string,"_vp.bin");
        vp_files[i] = open_bin_read(tmp_string);

        /* open vs file */ 
        strcpy(tmp_string,fname);
        strcat(tmp_string,"_vs.bin");
        vs_files[i] = open_bin_read(tmp_string);
        
        /* open rho file */ 
        strcpy(tmp_string,fname);
        strcat(tmp_string,"_rho.bin");
        rho_files[i] = open_bin_read(tmp_string);

        /* get number of values in each file */
        points_per_proc[i] = get_num_records(x_files[i]);
    } 

    /* allocate ibool model */
    allocate_ibool_model(imod,NPROC,points_per_proc);

    /* read values from bins */ 
    int points_read = 0;
    for (i=0;i<NPROC;i++)
    {
        read_bin_float(x_files[i]   , (imod->x)   + points_read);
        read_bin_float(z_files[i]   , (imod->z)   + points_read);
        read_bin_float(vp_files[i]  , (imod->vp)  + points_read);
        read_bin_float(vs_files[i]  , (imod->vs)  + points_read);
        read_bin_float(rho_files[i] , (imod->rho) + points_read);

        points_read += points_per_proc[i];
    }


    /* free local arrays */
    free(points_per_proc);


    /* close all files */
    for (i=0;i<NPROC;i++)
    {
        fclose(x_files[i]);
        fclose(z_files[i]);
        fclose(ibool_files[i]);
        fclose(vp_files[i]);
        fclose(vs_files[i]);
        fclose(rho_files[i]);
    }


    /* Free file pointer arrays */
    free(x_files);
    free(z_files);
    free(ibool_files);
    free(vp_files);
    free(vs_files);
    free(rho_files);
}

/*****************************************************************************/

void write_ibool_model_vpvs(ibool_model *imod, char *path)
{
    if (!imod->vpvs)
    {
        /* convert to vp/vs */
        fprintf(stderr,"Model is rmk\n");
        fprintf(stderr,"Conversion not yet implemented for ibool models\n");
        exit(1);
    }
    if (imod->natlog)
    {
        /* convert from natlog to regular */
        fprintf(stderr,"Model is natural log\n");
        fprintf(stderr,"Conversion not yet implemented for ibool models\n");
        exit(1);
    }


    int i;
    char fname[MAX_FNAME];
    char tmp_string[MAX_FNAME];
    int NPROC = imod->NPROC;


    FILE *x_file;
    FILE *z_file;
    FILE *vp_file;
    FILE *vs_file;
    FILE *rho_file;

    /* open files for writing and get array sizes*/
    int pts_written = 0;
    for (i=0;i<NPROC;i++)
    {
        /* write processor number to string */
        char proc_no[12];
        sprintf(proc_no,"%06d",i);

        /* write base file name */
        strcpy(fname,path);
        strcat(fname,"/proc");
        strcat(fname,proc_no);

        /* open x file */ 
        strcpy(tmp_string,fname);
        strcat(tmp_string,"_x.bin");
        x_file = open_bin_write(tmp_string);
        write_bin_float(imod->points[i], imod->x + pts_written, x_file);
        fclose(x_file);

        /* open z file */
        strcpy(tmp_string,fname);
        strcat(tmp_string,"_z.bin");
        z_file = open_bin_write(tmp_string);
        write_bin_float(imod->points[i], imod->z + pts_written, z_file);
        fclose(z_file);

        /* open vp file */ 
        strcpy(tmp_string,fname);
        strcat(tmp_string,"_vp.bin");
        vp_file = open_bin_write(tmp_string);
        write_bin_float(imod->points[i], imod->vp + pts_written, vp_file);
        fclose(vp_file);

        /* open vs file */ 
        strcpy(tmp_string,fname);
        strcat(tmp_string,"_vs.bin");
        vs_file = open_bin_write(tmp_string);
        write_bin_float(imod->points[i], imod->vs + pts_written, vs_file);
        fclose(vs_file);
        
        /* open rho file */ 
        strcpy(tmp_string,fname);
        strcat(tmp_string,"_rho.bin");
        rho_file = open_bin_write(tmp_string);
        write_bin_float(imod->points[i], imod->rho + pts_written, rho_file);
        fclose(rho_file);

        pts_written += imod->points[i];

    } 
}




