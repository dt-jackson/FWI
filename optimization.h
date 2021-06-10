#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H


#define WOLFE_C1 1e-4
#define WOLFE_C2 0.9
#define ZOOM_IT 20
#define ZOOM_TOL 1e-14


/*****************************************************************************/
/* typedefs, enums, and structs                                              */


typedef enum search_dir
{
    STEEPEST = 0,   /* steepest descent (negative gradient) direction */
    CG,             /* conjugate gradient direction (FR-PR method)    */
    LBFGS,          /* limited memory BFGS direction (quasi-Newton)   */
    NEWTON,         /* Newton direction - requires Hessian            */
} search_dir;


typedef struct optimization_parameters
{
    /* general */
    /* optimization stops when max_iterations have been performed or        */
    /* the gradient norm is less than stopping_tolerance                    */
    int max_iterations;
    double stopping_tolerance;


    /* search direction */
    search_dir sd;        /* 0=steepest, 1=conjugate grad, 2=LBFGS          */ 
    int LBFGS_mem;        /* number of models/gradients to store for L-BFGS */
    int angle_restart_allowed;  /* 1 if angle restart allowed, 0 if not     */


    /* step length selection algorithm */
    int backtrack_only;     /* 1 to perform backtracking only, 0 otherwise  */
    int bracket_only;       /* 1 to perform bracketing only, 0 otherwise    */
    int bracket_zoom_iterations;
    double bracket_zoom_tolerance;
    double initial_step_length;
    
    /* Wolfe conditions */ 
    int strong_curvature;   /* 1 to check strong curvature, 0 for regular   */
    double wolfe_c1;        /* parameter for sufficient decrease. Typ. 1e-4 */
    double wolfe_c2;        /* parameter for curvature. Typ. 0.5 - 0.9      */

    /* functions */
    double (*objective)(const double *x,const int N, void *aux);               /* function to calculate misfit   */
    void (*grad)(const double *x, double *gradient,const int N, void *aux);    /* function to calculate gradient */
    void (*hess)(const double *x, double *hessian, const int N, void *aux);    /* function to calculate hessian  */


} optimization_parameters;



/*****************************************************************************/
/* these functions perform the whole optimization procedure                  */
/* these are the ones that should be called by the user                      */



void line_search(double *x_in, 
                 double *x_out,
                 const int N,
                 optimization_parameters *opt_par,
                 void *aux
                 );




/*****************************************************************************/
/* dot products and norm                                                     */

double vec_dot_prod(const double *v1, const double *v2, const int N);
double vec_norm(const double *v1, const int N);


/*****************************************************************************/
/* estimate an appropriate initial step length for line search               */
/* not for the first iteration though                                        */

double initial_step_length(search_dir sd,
                           double alpha_old,
                           double dir_deriv_old,
                           double dir_deriv_new
                           );


/*****************************************************************************/
/* check if search direction is close to orthogonal to gradient              */

int check_angle_restart(double *grad, double *p, int N);




/*****************************************************************************/
/* search directions                                                         */

void steepest_descent_direction(const double *grad, double *p, const int N);

void newton_direction(const double *gradient, const double *hessian, const int N, double *p_new);


void conjugate_grad_direction(double *grad_old,
                              double *grad_new,
                              double *p_old,
                              double *p_new,
                              int N
                              );

void LBFGS_direction(double *grad_new,
                     double *grad_old,
                     double *model_new,
                     double *model_old,
                     double *p_new,
                     int N,
                     int m,
                     int k
                     );

/*****************************************************************************/

#define HESSMOD_MAX_IT 100
int hessian_modification(const double *hessian, double *upper_factor, int N, double beta);
void modified_chol(double *A, int N);
/*****************************************************************************/
/* Wolfe conditions                                                          */
/* return 1 if condition is satsified, 0 otherwise                           */

int check_armijo(double new_misfit, double old_misfit, 
                 double dir_deriv, double step_len);
int check_curvature(double dir_deriv_old, double dir_deriv_new);
int check_strong_curvature(double dir_deriv_old, double dir_deriv_new);



/*****************************************************************************/
/* bracketing line search functions*/

double ls_cubic_interpolation(const double *x, 
               const double *p,
               const int N,
               optimization_parameters *opt_par,
               const double bound1,
               const double bound2,
               void *aux
               );

double ls_zoom(const double *x, 
               const double *p,
               const int N,
               optimization_parameters *opt_par,
               double a_lo,
               double a_hi,
               void *aux
               );

double ls_bracket(const double *x, 
               const double *p,
               const int N,
               optimization_parameters *opt_par,
               double alpha0,
               void *aux
               );


/*****************************************************************************/
/* backtracking line search functions */

double back_track_quad(double phi0, double phip0, double phi_a0, double a0);

double back_track_cubic(double phi0, 
                        double phip0,
                        double phi_a0,
                        double phi_a1,
                        double a0,
                        double a1
                        );

double back_track(double *x,
                  double *p,
                  int N,
                  optimization_parameters *opt_par,
                  double alpha,
                  void *aux
                  );

void* opt_allocate(long nbytes);

#endif