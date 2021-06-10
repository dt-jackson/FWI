#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "read_parameters.h"



par_struct* setup_parameters(parameters *p, par_struct *pp)
{
    /* to add parameters, add there in the parameters structure, 
    update the function setup_parameters, and adjust the num_par macro  */


    
    int ind;

    /* DT */
    ind = 0; 
    strcpy(pp[ind].par_name,"DT");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->dt;

    /* NSTEP */
    ind++;
    strcpy(pp[ind].par_name,"NSTEP");
    pp[ind].pt = INTEGER;
    pp[ind].par = &p->nstep;

    /* NUM_SRC */   
    ind++;
    strcpy(pp[ind].par_name,"NUM_SRC");
    pp[ind].pt = INTEGER;
    pp[ind].par = &p->num_src;

    /* NUM_REC */   
    ind++;
    strcpy(pp[ind].par_name,"NUM_REC");
    pp[ind].pt = INTEGER;
    pp[ind].par = &p->num_rec;
    
    /* NPROC */   
    ind++;
    strcpy(pp[ind].par_name,"NPROC");
    pp[ind].pt = INTEGER;
    pp[ind].par = &p->nproc;

    /* XMIN */   
    ind++;
    strcpy(pp[ind].par_name,"XMIN");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->xmin;

    /* XMAX */   
    ind++;
    strcpy(pp[ind].par_name,"XMAX");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->xmax;
   
    /* ZMIN */   
    ind++;
    strcpy(pp[ind].par_name,"ZMIN");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->zmin;

    /* ZMAX */   
    ind++;
    strcpy(pp[ind].par_name,"ZMAX");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->zmax;
    
    /* DX_REG */   
    ind++;
    strcpy(pp[ind].par_name,"DX_REG");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->dx_reg;
   
    /* DX_INV */   
    ind++;
    strcpy(pp[ind].par_name,"DX_INV");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->dx_inv;

    /* SP_LP_CUTOFF */   
    ind++;
    strcpy(pp[ind].par_name,"SP_LP_CUTOFF");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->sp_lp_cutoff;

    /* MED_FILT_SIZE */   
    ind++;
    strcpy(pp[ind].par_name,"MED_FILT_SIZE");
    pp[ind].pt = INTEGER;
    pp[ind].par = &p->med_filt_size;

    /* SEARCH_DIR */   
    ind++;
    strcpy(pp[ind].par_name,"SEARCH_DIR");
    pp[ind].pt = INTEGER;
    pp[ind].par = &p->search_direction;

    /* INITIAL_STEP_MULT */   
    ind++;
    strcpy(pp[ind].par_name,"INITIAL_STEP_MULT");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->init_step_mult;

    /* MODEL_PATH */   
    ind++;
    strcpy(pp[ind].par_name,"MODEL_PATH");
    pp[ind].pt = STRING;
    pp[ind].par = p->model_path;

    /* SF_BIN_FORWARD */   
    ind++;
    strcpy(pp[ind].par_name,"SF_BIN_FORWARD");
    pp[ind].pt = STRING;
    pp[ind].par = p->sf_bin_forward;

    /* SF_BIN_ADJOINT */   
    ind++;
    strcpy(pp[ind].par_name,"SF_BIN_ADJOINT");
    pp[ind].pt = STRING;
    pp[ind].par = p->sf_bin_adjoint;

    /* SRC_LOC_FILE */   
    ind++;
    strcpy(pp[ind].par_name,"SRC_LOC_FILE");
    pp[ind].pt = STRING;
    pp[ind].par = p->src_loc_file;

    /* REC_LOC_FILE */   
    ind++;
    strcpy(pp[ind].par_name,"REC_LOC_FILE");
    pp[ind].pt = STRING;
    pp[ind].par = p->rec_loc_file;

    /* MIN_SR_DIST */   
    ind++;
    strcpy(pp[ind].par_name,"MIN_SR_DIST");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->min_sr_dist;

    /* MAX_SR_DIST */   
    ind++;
    strcpy(pp[ind].par_name,"MAX_SR_DIST");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->max_sr_dist;
    
    
    
    /* VP_MAX */   
    ind++;
    strcpy(pp[ind].par_name,"VP_MAX");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->max_vp;

    /* VP_MIN */   
    ind++;
    strcpy(pp[ind].par_name,"VP_MIN");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->min_vp;

    /* VS_MAX */   
    ind++;
    strcpy(pp[ind].par_name,"VS_MAX");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->max_vs;

    /* VS_MIN */   
    ind++;
    strcpy(pp[ind].par_name,"VS_MIN");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->min_vs;


    /* RHO_MAX */   
    ind++;
    strcpy(pp[ind].par_name,"RHO_MAX");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->max_rho;

    /* RHO_MIN */   
    ind++;
    strcpy(pp[ind].par_name,"RHO_MIN");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->min_rho;

    /* PR_MAX */   
    ind++;
    strcpy(pp[ind].par_name,"PR_MAX");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->max_pr;

    /* PR_MIN */   
    ind++;
    strcpy(pp[ind].par_name,"PR_MIN");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->min_pr;

    /* MASK_DENS */   
    ind++;
    strcpy(pp[ind].par_name,"MASK_DENS");
    pp[ind].pt = INTEGER;
    pp[ind].par = &p->mask_dens;

    /* MASK_KERNELS */   
    ind++;
    strcpy(pp[ind].par_name,"MASK_KERNELS");
    pp[ind].pt = INTEGER;
    pp[ind].par = &p->mask_kernels;


    /* MASK_XMIN */   
    ind++;
    strcpy(pp[ind].par_name,"MASK_XMIN");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->mask_xmin;
    
    /* MASK_XMAX */   
    ind++;
    strcpy(pp[ind].par_name,"MASK_XMAX");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->mask_xmax;
    
    /* MASK_ZMIN */   
    ind++;
    strcpy(pp[ind].par_name,"MASK_ZMIN");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->mask_zmin;
    
    /* MASK_ZMAX */   
    ind++;
    strcpy(pp[ind].par_name,"MASK_ZMAX");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->mask_zmax;



    /* MAX_ITERATIONS */   
    ind++;
    strcpy(pp[ind].par_name,"MAX_ITERATIONS");
    pp[ind].pt = INTEGER;
    pp[ind].par = &p->max_iterations;


    /* STOPPING_TOL */   
    ind++;
    strcpy(pp[ind].par_name,"STOPPING_TOL");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->stopping_tol;

    /* LBFGS_MEM */   
    ind++;
    strcpy(pp[ind].par_name,"LBFGS_MEM");
    pp[ind].pt = INTEGER;
    pp[ind].par = &p->LBFGS_mem;

    /* ANGLE_RESTART */   
    ind++;
    strcpy(pp[ind].par_name,"ANGLE_RESTART");
    pp[ind].pt = INTEGER;
    pp[ind].par = &p->angle_restart;

    /* BACKTRACK_ONLY */   
    ind++;
    strcpy(pp[ind].par_name,"BACKTRACK_ONLY");
    pp[ind].pt = INTEGER;
    pp[ind].par = &p->backtrack_only;

    /* BRACKET_ONLY */   
    ind++;
    strcpy(pp[ind].par_name,"BRACKET_ONLY");
    pp[ind].pt = INTEGER;
    pp[ind].par = &p->bracket_only;

    /* STRONG_CURVATURE */   
    ind++;
    strcpy(pp[ind].par_name,"STRONG_CURVATURE");
    pp[ind].pt = INTEGER;
    pp[ind].par = &p->strong_curvature;

    /* WOLFE_C1 */   
    ind++;
    strcpy(pp[ind].par_name,"WOLFE_C1");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->wolfe_c1;

    /* WOLFE_C2 */   
    ind++;
    strcpy(pp[ind].par_name,"WOLFE_C2");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->wolfe_c2;


    /* DATA_STDEV */   
    ind++;
    strcpy(pp[ind].par_name,"DATA_STDEV");
    pp[ind].pt = DOUBLE;
    pp[ind].par = &p->data_stdev;

    
    /* MISFIT_TYPE */   
    ind++;
    strcpy(pp[ind].par_name,"MISFIT_TYPE");
    pp[ind].pt = INTEGER;
    pp[ind].par = &p->misfit_type;



    p->num_parameters = ind+1;
    if (ind+1 > NUM_PAR)
    {
        fprintf(stderr,"NUM_PAR too small. Increase and recompile\n");
        exit(1);
    }

    return pp;

    /* to add parameters, add there in the parameters structure, 
    update the function setup_parameters, and adjust the num_par macro  */
}




void read_parameters(parameters *p, char *path)
{
    par_struct *pp = malloc(sizeof(*pp) * NUM_PAR);
    if (!pp)
    {
        fprintf(stderr,"Error: setup_parameters: memory allocation failure\n");
        exit(1);
    }
    
    pp = setup_parameters(p, pp);



    FILE *fid = fopen(path,"r");
    if (!fid)
    {
        fprintf(stderr,"Error reading file %s\n",path);
        exit(1);
    }
    char *line;
    size_t line_len = 0;

    char *par_name = NULL;
    char *par_value = NULL;

    int i;
    /* read each line and look for parameter names */ 
    while (getline(&line,&line_len,fid) != -1)
    {
        
        /* remove comments */
        remove_comments(line);
        
        ///* remove white space */
        remove_white_space(line);


        if (get_tokens(line, &par_name, &par_value))
        {
            /* make sure the parameter name is upper case */            
            line_to_upper(par_name);

            /* check if the par_name matches */
            for (i=0;i<p->num_parameters;i++)
            {
                if (strcmp(par_name, pp[i].par_name) == 0)
                {
                    /* if there is a match, look at the type and
                    store it in the appropriate location */
                    switch (pp[i].pt)
                    {
                        case INTEGER:
                            *(int*)pp[i].par = atoi(par_value);
                            break;
                        case DOUBLE:
                            *(double*)pp[i].par = atof(par_value);
                            break;
                        case STRING:
                            strcpy((char*)pp[i].par,par_value);
                            break;
                        default:
                            fprintf(stderr,"Parameter type not found\n");
                            exit(1);
                    }
                    break;
                }
            }

        }

    }
    
    free(line);
    free(pp);
    fclose(fid);
}




int get_tokens(char *line, char **par_name, char **par_value)
{
    char delimiter[] = "=";
    *par_name = strtok(line,delimiter);
    *par_value = strtok(NULL,delimiter);

    if (*par_value == NULL)
    {
        /* no delimiter found */
        return 0;
    }
    return 1;
}




void remove_comments(char *line)
{
    /* will remove bash style comments (begin with # going to end of line) */
    char *comment = NULL;
    comment = strchr(line, '#');
    if (comment)
        *comment = '\0';
}


void remove_white_space(char *line)
{
    char temp[LINE_MAX];
    int i;
    int len = strlen(line);

    int count = 0;
    for (i=0;i<len;i++)
    {
        if(line[i] != ' ' && line[i] != '\t' && line[i] != '\n')
            temp[count++] = line[i];
    }
    temp[count] = '\0';

    for (i=0;i<count+1;i++)
    {
        line[i] = temp[i];
    }
}

void line_to_upper(char *line)
{
    int len = strlen(line);
    int i;
    for (i=0;i<len;i++)
        line[i]  = toupper(line[i]);
}




void print_parameters(parameters *p)
{
    printf("dt: %f\n",p->dt);
    printf("nstep: %d\n",p->nstep);
    printf("nproc: %d\n",p->nproc);

    printf("data_stdev: %E\n",p->data_stdev);

    printf("num_src: %d\n",p->num_src);
    printf("num_rec: %d\n",p->num_rec);
    printf("Source location file: %s\n",p->src_loc_file);
    printf("Receiver location file: %s\n",p->rec_loc_file);
    printf("Min S-R distance: %f\n", p->min_sr_dist);
    printf("Max S-R distance: %f\n", p->max_sr_dist);


    printf("xmin: %f\n",p->xmin);
    printf("xmax: %f\n",p->xmax);
    printf("zmin: %f\n",p->zmin);
    printf("zmin: %f\n",p->zmax);

    printf("reg grid size: %f\n",p->dx_reg);
    printf("inv grid size: %f\n",p->dx_inv);

    printf("sptial filter cutoff: %f\n",p->sp_lp_cutoff);
    printf("median filter size: %d\n",p->med_filt_size);

    printf("search direction: %d\n",p->search_direction);
    printf("init step len: %f\n",p->init_step_mult);

    printf("model path: %s\n",p->model_path);
    printf("sf forward path: %s\n",p->sf_bin_forward);
    printf("sf adjoint path: %s\n",p->sf_bin_adjoint);
    
    printf("mask density: %d\n", p->mask_dens);
    printf("mask kernels: %d\n", p->mask_kernels);

    printf("mask xmin: %f\n", p->mask_xmin);
    printf("mask xmax: %f\n", p->mask_xmax);
    printf("mask zmin: %f\n", p->mask_zmin);
    printf("mask zmax: %f\n", p->mask_zmax);

    printf("Vp min: %f\n", p->min_vp);
    printf("Vp max: %f\n", p->max_vp);
    printf("Vs min: %f\n", p->min_vs);
    printf("Vs max: %f\n", p->max_vs);
    printf("Density min: %f\n", p->min_rho);
    printf("Density max: %f\n", p->max_rho);
    printf("Poisson ratio min: %f\n", p->min_pr);
    printf("Poisson ratio max: %f\n", p->max_pr);
}







