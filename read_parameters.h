#ifndef READ_PARAMETERS_H
#define READ_PARAMETERS_H

#ifndef LINE_MAX
#define LINE_MAX 512
#endif

#ifndef PARNAME_MAX
#define PARNAME_MAX 50
#endif

#ifndef FILENAME_MAX
#define FILENAME_MAX 256
#endif

#define NUM_PAR 150 /* <- this is the max. If there are more than 150 parameters, increase this */



typedef enum par_type
{
    INTEGER = 0,
    DOUBLE = 1,
    STRING = 2,
} par_type;



typedef struct par_struct
{
    char par_name[PARNAME_MAX];
    par_type pt;
    void *par;
} par_struct;




typedef struct parameters
{
    /* To add parameters, add there in the parameters structure, 
    update the function setup_parameters, and adjust the NUM_PAR macro  */

    /* simulation parameters */
    double dt;
    int nstep;

    /* sources, receivers, processors */
    int num_src;
    int num_rec;
    int nproc;

    /* grid extent */
    double xmin;
    double xmax;
    double zmin;
    double zmax;

    /* grid sizes (same in x and z directions) */
    double dx_reg;
    double dx_inv;

    /* spatial lowpass filter cutoff frequency */
    double sp_lp_cutoff;
    
    /* Tikhonov regularization parameters */ 
    //double tr_lambda_x;
    //double tr_lambda_z;

    /* optimization parameters */
    int search_direction;
    double init_step_mult;

    /* directory where model files are located */
    char model_path[FILENAME_MAX];

    char sf_bin_forward[FILENAME_MAX];
    char sf_bin_adjoint[FILENAME_MAX];

    /* filenames for src and rec locations */
    char src_loc_file[FILENAME_MAX];
    char rec_loc_file[FILENAME_MAX];

    double min_sr_dist;
    double max_sr_dist;
    

    double max_vp;
    double min_vp;
    double max_vs;
    double min_vs;
    double max_rho;
    double min_rho;
    double min_pr;
    double max_pr;

    int med_filt_size;


    int mask_dens;

    int mask_kernels;
    double mask_xmin;
    double mask_xmax;
    double mask_zmin;
    double mask_zmax;


    int max_iterations;
    double stopping_tol;
    int LBFGS_mem;
    int angle_restart;
    int backtrack_only;
    int bracket_only;
    int strong_curvature;
    double wolfe_c1;
    double wolfe_c2;
    int misfit_type;

    double data_stdev;


    int num_parameters; /* total number of parameters */

    /* To add parameters, add there in the parameters structure, 
    update the function setup_parameters */


} parameters;

void read_parameters(parameters *p, char *path);

void remove_comments(char *line);

void remove_white_space(char *line);

int get_tokens(char *line, char **par_name, char **par_value);

par_struct* setup_parameters(parameters *p, par_struct *pp);

void print_parameters(parameters *p);

void line_to_upper(char *line);


#endif