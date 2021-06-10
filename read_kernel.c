#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "unique_model.h"
#include "ibool_model.h"
#include "bin_io.h"

#define MAX_FILE_LEN 256


void read_event_kernel_rmk(ibool_model *imod, char *path, int NPROC)
{
    FILE *rho_file;
    FILE *mu_file;
    FILE *kappa_file;

    imod->vpvs = 0;
    imod->kernel = 1;
    imod->natlog = 0; /* not really valid for kernels */

    char proc_no[12];
    char fname[MAX_FILE_LEN];
    fname[0] = '\0';

    /* read rho kernel */
    int ii, jj;
    int points_read = 0;
    for (ii=0; ii<NPROC; ii++)
    {
        sprintf(proc_no,"%06d",ii);
        strcat(fname,path);
        strcat(fname,"/proc");
        strcat(fname,proc_no);
        strcat(fname,"_rho_kernel.bin");
        rho_file = open_bin_read(fname);
        read_bin_float(rho_file,imod->rho + points_read);
        
        points_read += imod->points[ii];
        fclose(rho_file);
        fname[0] = '\0';
    }
    
    
    /* read mu kernel */
    points_read = 0;
    for (ii=0; ii<NPROC; ii++)
    {
        sprintf(proc_no,"%06d",ii);
        strcat(fname,path);
        strcat(fname,"/proc");
        strcat(fname,proc_no);
        strcat(fname,"_mu_kernel.bin");
        mu_file = open_bin_read(fname);
        read_bin_float(mu_file,imod->mu + points_read);
        
        points_read += imod->points[ii];
        fclose(mu_file);
        fname[0] = '\0';
    }
    
    
    /* read kappa kernel */
    points_read = 0;
    for (ii=0; ii<NPROC; ii++)
    {
        sprintf(proc_no,"%06d",ii);
        strcat(fname,path);
        strcat(fname,"/proc");
        strcat(fname,proc_no);
        strcat(fname,"_kappa_kernel.bin");
        kappa_file = open_bin_read(fname);
        read_bin_float(kappa_file,imod->kappa + points_read);
        
        points_read += imod->points[ii];
        fclose(kappa_file);
        fname[0] = '\0';
    }


}




