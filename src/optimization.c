#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>

#include "optimization.h"







/*****************************************************************************/
/* vector routines */

double vec_dot_prod(const double *v1, const double *v2, const int N)
{
    int i;
    double nn = 0;
    for (i=0;i<N;i++)
    {
        nn += v1[i] * v2[i];
    }
    return nn;
}

double vec_norm(const double *v1, const int N)
{
    return sqrt(vec_dot_prod(v1, v1, N));
}


/*****************************************************************************/
/* Wolfe conditions */

int check_armijo(double new_misfit, double old_misfit, 
                 double dir_deriv, double step_len)
{
    /* returns 1 if Armijo condition satisfied, 0 otherwise */
    if (new_misfit <= old_misfit + WOLFE_C1 * step_len * dir_deriv) 
    {
        /* step length okay */
        return 1;
    }
    return 0;
}

int check_curvature(double dir_deriv_old, double dir_deriv_new)
{
    /* returns 1 if curvature satisfied, 0 otherwise */
    if (dir_deriv_new >= WOLFE_C2 * dir_deriv_old)
    {
        return 1;
    }
    return 0;
}

int check_strong_curvature(double dir_deriv_old, double dir_deriv_new)
{
    /* returns 1 if strong curvature satisfied, 0 otherwise */
    if (dir_deriv_new >=  WOLFE_C2 * dir_deriv_old && 
        dir_deriv_new <= -WOLFE_C2 * dir_deriv_old)
    {
        return 1;
    }
    return 0;
}


/*****************************************************************************/

double initial_step_length(search_dir sd, double alpha_old, double dir_deriv_old, double dir_deriv_new)
{
    /* calculate initial step length for a line search step         */
    /* inputs:                                                      */
    /*  sd: search direction, CG, STEEPEST, or LBFGS                */
    /*  alpha_old: previous (accepted) step length                  */
    /*  dir_deriv_old: directional derivative at previous step      */
    /*      (dot product of previous gradient and search direction) */
    /*  dir_deriv_new: directional derivative at current step       */
    /*      (dot product of current gradient and search direction)  */
    


    if (sd == CG || sd == STEEPEST)
    {
        /* see Nocedal  & Wright 2006, p. 59 */
        return alpha_old * (dir_deriv_old / dir_deriv_new);
    }

    /* quasi-Newton methods are well scaled use steplength of 1.0 for L-BFGS */
    return 1.0;

}



/*****************************************************************************/
/* steepest descent direction */

void steepest_descent_direction(const double *grad, double *p, const int N)
{
    int i;
    for (i=0;i<N;i++)
    {
        p[i] = -grad[i];
    }
}

/*****************************************************************************/
/* conjugate gradient direction */

void conjugate_grad_direction(double *grad_old, double *grad_new, double *p_old, double *p_new, int N)
{
    /* 
    Determine a search direction by the conjugate gradient method

    Uses PR-FR conjugate gradient method to calculate a conjugate direection,
    see Nocedal & Wright 2006, p. 123. 

    All vectors should be double precision of length N.
    grad_old: previous gradient
    grad_new: current gradient
    p_old: previous search direction
    p_new: new search direction (output, must be allocated already)
    */

    /* calculate beta values for both PR and FR methods */
    double beta_lower = 0;
    double beta_upper_FR = 0;
    double beta_upper_PR = 0;

    int i;
    for (i=0;i<N;i++)
    {
        beta_lower += grad_old[i] * grad_old[i];
        beta_upper_FR += grad_new[i] * grad_new[i];
        beta_upper_PR += grad_new[i] * (grad_new[i] - grad_old[i]);
    }

    double beta_FR = beta_upper_FR / beta_lower;
    double beta_PR = beta_upper_PR / beta_lower;


    /* select beta based on criteria given in Nocedal and Wright 2006 */
    double beta = 0;

    if (beta_PR < -beta_FR)
    {
        beta = -beta_FR;
    }
    else if (fabs(beta_PR) <= beta_FR)
    {
        beta = beta_PR;
    }
    else if (beta_PR > beta_FR)
    {
        beta = beta_FR;
    }

    /* calculate new search direction */
    for (i=0;i<N;i++)
    {
        p_new[i] = -grad_new[i] + beta * p_old[i];
    }
}



/*****************************************************************************/
/* L-BFGS direction */


void LBFGS_direction(double *grad_new, double *grad_old, double *model_new, double *model_old, double *p_new, int N, int m, int k)
{
    /* call with the last parameter negative to clean up local variables */

    int ii, jj;

    static double *s = NULL; 
    static double *y = NULL;
   
    static double *q = NULL;
    static double *aa = NULL;
    static double *rho = NULL;
    static int *ind = NULL;

    double *r = p_new;


    /* allocate memory */
    /* note: opt_allocate is the same as malloc, but checks to make sure    */ 
    /* memory allocate was successful                                       */
    if (!s)
    {
        s = opt_allocate(sizeof(*s) * N * m);
    }

    if (!y)
    {
        y = opt_allocate(sizeof(*y) * N * m);
    }

    if (!q)
    {
        q = opt_allocate(sizeof(*q) * N);
    }

    if (!aa)
    {
        aa = opt_allocate(sizeof(*aa) * m);
    }

    if (!rho)
    {
        rho = opt_allocate(sizeof(*rho) * m);
    }

    if (!ind)
    {
        ind = opt_allocate(sizeof(*ind) * m);
        for (ii=0;ii<m;ii++)
        {
            ind[ii] = ii * N;
        }    
    }

    /* j is the min of k and m */
    int j = (k < m) ? k : m;

    /* position to write new y and s vectors */ 
    int pos = ((k-1) % m) * N;

    /* write new y and s vectors */
    for (ii=0;ii<N;ii++)
    {
        y[ii+pos] = grad_new[ii] - grad_old[ii];
        s[ii+pos] = model_new[ii] - model_old[ii];
    }

    /* calculate rho */
    rho[pos / N] = 1.0 / vec_dot_prod(y+pos, s+pos, N);


    /* calculate gamma */
    double gamma = vec_dot_prod(s+pos, y+pos, N) / vec_dot_prod(y+pos, y+pos, N);

    if (k<m)
    {
        /* shift ind vector so ind[0] is still the index 
        of the oldest element */
        int first_ele = ind[0];
        for (ii=0;ii<m-1;ii++)
        {
            ind[ii] = ind[ii+1];
        }
        ind[m-1] = first_ele;
    }
    
    
    /* assign initial q value */
    for (jj=0;jj<N;jj++)
    {
        q[jj] = grad_new[jj];
    }

    /*begin backwards loop */
    for (ii=j-1;ii>=0;ii--)
    {
        aa[ii] = rho[ind[ii] / N] * vec_dot_prod(s+ind[ii], q, N);


        for (jj=0;jj<N;jj++)
        {
            q[jj] = q[jj] - aa[ii] * (y+ind[ii])[jj];
        }
    }

    /* assign r */
    for (jj=0;jj<N;jj++)
    {
        r[jj] = gamma * q[jj];
    }


    /* begin forwards loop */
    double beta;
    for (ii=0;ii<j;ii++)
    {
        beta = rho[ind[ii] / N] * vec_dot_prod(y+ind[ii], r, N);


        for (jj=0;jj<N;jj++)
        {
            r[jj] = r[jj] + (s+ind[ii])[jj] * (aa[ii] - beta);
        }
    }

    /* negate r */
    for (jj=0;jj<N;jj++)
    {
        r[jj] = -r[jj];
    }

}

/*****************************************************************************/
#include <lapacke.h>
void newton_direction(const double *gradient, const double *hessian, const int N, double *p_new)
{
    /********************************************************/
    /* not implemented yet - use steepest descent for now   */
    //fprintf(stderr,"Newton direction not implemented yet\n");
    //fprintf(stderr,"Reverting to steepest descent\n");
    //steepest_descent_direction(gradient, p_new, N);
    /********************************************************/

    /* make a temp copy of the gradient*/
    double *temp_g = opt_allocate(sizeof(*temp_g) * N);
    double *temp_h = opt_allocate(sizeof(*temp_h) * N * N);

    memcpy(temp_g, gradient, sizeof(*temp_h) * N );


    /* need to check if Hessian is PD and modify if not                 */
    /* Note that this repeated performs a Cholesky factorization on     */
    /* the Hessian and will be expensive if N is large                  */
    /* In those cases, a modified Cholesky factorization (see Nocedal   */
    /* and Wright 2006, pp. 52-54) may be a better approach             */
    double beta = sqrt(DBL_EPSILON);
    int pd = hessian_modification(hessian, temp_h, N, beta);

    /* if we found a pd Hessian approximation, solve system for p       */
    /* otherwise, return steepest descent direction                     */
    if (pd != 1)
    {
        fprintf(stderr,"Could not find a positive definite "
                       "Hessian approximation\n");
        fprintf(stderr,"Returning steepest descent direction\n");
    }
    else
    {
        /* solve system for p: Hp = -g */
        int nrhs = 1;
        int lda = N;
        int ldb = N;
        LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', N, nrhs, temp_h, lda, temp_g, ldb);
    }

    /* negate gradient or solution to system */
    int i;
    for (i = 0;i<N;i++)
        p_new[i]  = -temp_g[i];
    
    free(temp_h);
    free(temp_g);
}

#define HESSMOD_MAX_IT 100
int hessian_modification(const double *hessian, double *upper_factor, int N, double beta)
{
    /* beta may need to be adjusted based on the scaling of the problem */

    /* find minimum diagonal element of the hessian */
    int i;
    double min_diag = hessian[0];
    for (i=0;i<N;i++)
    {
        double diag = hessian[i*(N+1)];
        min_diag = (diag < min_diag) ? diag : min_diag;
    }

    /* find initial tau */
    double tau;
    if (min_diag > 0)
    {
        /* if it is already pd, min_diag will be positive */
        tau = 0.0;
    }
    else 
    {
        tau = - min_diag + beta;
    }

    /* add multiples of the identity until hessian is positive definite */
    int k;
    int pd;
    int lda = N;
    for (k=0;k<HESSMOD_MAX_IT;k++)
    {
        /* copy hessian and add identity multiple */
        /* upper_factor will change every loop because of the dpotrf call */
        memcpy(upper_factor, hessian, sizeof(*hessian) * N *N);
        for (i=0;i<N;i++)
        {
            upper_factor[i * (N+1)] = hessian[i * (N+1)] + tau;
        } 

        /* try cholesky factorization */
        pd = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', N, upper_factor, lda);
        if (pd == 0)
        {
            /* positive definite */
            return 1;
        }
        /* increment tau */
        tau = 2.0 * tau > beta ? 2.0 * tau : beta;

    }
    /* return 0 if not positive definite */
    return 0;
    
}


/*****************************************************************************/
/* function for performing line search */


void line_search(double *x_in, 
                       double *x_out,
                       const int N,
                       optimization_parameters *opt_par,
                       void *aux
                       )
{

    double (*misfit)(const double *, const int, void *) = opt_par->objective;
    void (*grad)(const double *, double *, int, void *) = opt_par->grad;
    void (*hess)(const double *, double *, int, void *) = opt_par->hess;

    int success = 0;
    int i,j;

    /* note: opt_allocate is the same as malloc, but checks to make sure    */ 
    /* memory allocate was successful                                       */

    double *gradient = opt_allocate(sizeof(*gradient) * N);
    double *grad_old = opt_allocate(sizeof(*grad_old) * N);

    double *p = opt_allocate(sizeof(*p) * N);
    double *p_old = opt_allocate(sizeof(*p_old) * N);

    double *x_old = opt_allocate(sizeof(*x_old) * N);
    double *x_trial = opt_allocate(sizeof(*x_trial) * N);

    double *hessian = NULL;
    if (opt_par->sd == NEWTON)
        hessian = opt_allocate(sizeof(*hessian) * N * N);

    /* initial alpha_old and dir_deriv_old to keep the compiler from complaining,   */
    /* but these initialization values will never be used                           */
    double alpha_old = 1;
    double alpha_new;

    double dir_deriv_old = 1;
    double dir_deriv_new;

    /* copy x_in to x_old */
    for (i=0;i<N;i++)
    {
        x_out[i] = x_in[i];
    }


    /* evaluate function and gradient at initial point */
    double f_old = misfit(x_out, N, aux);
    grad(x_out, gradient, N, aux);

    /* check grad norm for exit condition */
    if (vec_norm(gradient, N) < opt_par->stopping_tolerance)
    {
        success = 1;
        free(p);
        free(gradient);
        return;
    }

    int m = opt_par->LBFGS_mem;

    /* start main loop */
    search_dir sd;
    int max_it = opt_par->max_iterations;
    for (j=0;j<max_it;j++)
    {
        printf("Starting line search, iteration %d\n",j);
        if (j == 0 && opt_par->sd != NEWTON)
        {
            /* first iteration - use steepest descent */
            /* get steepest descent direction */
            steepest_descent_direction(gradient, p, N);
            sd = STEEPEST;
            dir_deriv_new = vec_dot_prod(gradient, p, N);
            alpha_new = opt_par->initial_step_length;
        }
        else
        {
            if (opt_par->sd == LBFGS)
            {
                LBFGS_direction(gradient, grad_old, x_out, x_old, p, N, m, j);
                

                sd = LBFGS;
                dir_deriv_new = vec_dot_prod(gradient, p, N);
            }
            else if (opt_par->sd == CG)
            {
                conjugate_grad_direction(grad_old, gradient, p_old, p, N);
                sd = CG;
                dir_deriv_new = vec_dot_prod(gradient, p, N);
            }
            else if (opt_par->sd == NEWTON)
            {
                /* compute hessian */
                hess(x_out, hessian, N, aux);
                newton_direction(gradient, hessian, N, p);
                sd = NEWTON;
                dir_deriv_new = vec_dot_prod(gradient, p, N);
            }

            if (opt_par->sd == STEEPEST || (opt_par->angle_restart_allowed && check_angle_restart(gradient, p, N)))
            {
                steepest_descent_direction(gradient, p, N);
                sd = STEEPEST;
                dir_deriv_new = vec_dot_prod(gradient, p, N);
            }

            alpha_new = initial_step_length(sd, alpha_old, dir_deriv_old, dir_deriv_new);
        }

        if (opt_par->bracket_only || j==0)
        {
            alpha_new = ls_bracket(x_out, p, N, opt_par, alpha_new, aux);
        }
        else 
        {
            alpha_new = back_track(x_out, p, N, opt_par, alpha_new, aux);
        }
        

        for (i=0;i<N;i++)
        {
            p_old[i] = p[i];
            grad_old[i] = gradient[i];
            x_old[i] = x_out[i];
        }


        if (!opt_par->backtrack_only && !opt_par->bracket_only && j!=0)
        {
            /* update model */
            for (i=0;i<N;i++)
            {
                x_trial[i] = x_out[i] + alpha_new * p[i];
            }

            /* evaluate gradient at new point */
            grad(x_trial, gradient, N, aux);

            /* check (strong) curvature */ 
            double ppi = vec_dot_prod(gradient, p, N);
            if (!check_strong_curvature(dir_deriv_new, ppi))
            {
                //alpha_new = initial_step_length(sd, alpha_old, dir_deriv_old, dir_deriv_new);
                alpha_new = ls_bracket(x_out, p, N, opt_par, alpha_new, aux); 
            }

        }

        alpha_old = alpha_new;
        dir_deriv_old = dir_deriv_new;


        /* calculate new x value */
        for (i=0; i<N; i++)
        {
            x_out[i] = x_out[i] + alpha_new * p[i];
        }


        /* evaluate gradient at new point */
        grad(x_out, gradient, N, aux);


        /* check grad norm for exit condition */
        if (vec_norm(gradient, N) < opt_par->stopping_tolerance)
        {
            success = 1;
            break;
        }

    }

    
    printf("No. iterations: %d\n",j);

    free(p);
    free(p_old);
    free(gradient);
    free(grad_old);
    free(x_trial);
    free(x_old);

    if (hessian)
        free(hessian);

    if (!success)
    {
        fprintf(stderr,"Error: line search failed to find a solution\n");
        fprintf(stderr,"Point returned is not stationary\n");
    }

}







/*****************************************************************************/

double ls_zoom(const double *x, 
               const double *p,
               const int N,
               optimization_parameters *opt_par,
               double a_lo,
               double a_hi,
               void *aux
               )
{
    double (*misfit)(const double *, const int, void *) = opt_par->objective;
    void (*grad)(const double *, double *, int, void *) = opt_par->grad;


    int success = 0;

    double bb_star = 0.0;

    double mf;
    double *gradient = opt_allocate(sizeof(*gradient) * N);
    double *x_update = opt_allocate(sizeof(*x_update) * N);

    /* calculate initial misfit and gradient */
    mf = misfit(x, N, aux);
    grad(x, gradient, N, aux);


    double phi0 = mf;
    double phi_prime0 = vec_dot_prod(gradient, p, N);


    double phi_aj;
    double phi_prime_aj;

    double phi_a_lo;
    
    double mean_step;

    double tol = 0;

    /* start iterations */
    double a_j;
    int i, j;
    for (i=0;i<opt_par->bracket_zoom_iterations;i++)    
    {
        printf("Bracket successful, starting zoom, iteration %d\n",i);
        printf("a_lo, a_hi: %f, %f\n", a_lo, a_hi);
        fflush(stdout);

        /* check if a_lo and a_hi are too close together */
        mean_step = 0.5 * (a_lo + a_hi);
        if (fabs(a_lo - a_hi) / mean_step < 0.01)
        {
            bb_star = mean_step;
            success = 1;
            break;
        }
        if (mean_step < 1e-6) /* <-- this should really depend on the scaling of the problem */
        {
            success = 0;
            break;
        }

        a_j = ls_cubic_interpolation(x, p, N, opt_par, a_lo, a_hi, aux);
        //a_j = (a_lo + a_hi) / 2.0;

        tol = mean_step * 0.01; /* one percent of the mean step */
        //fprintf(stderr, "Trial step: %f\n",a_j);
        if (fabs(a_j - a_lo) < tol || fabs(a_j - a_hi) < tol)
        {
            /* too close to boundary of interval */
            //fprintf(stderr,"ls_zoom: too close to boundary\n");
            //success = 0;
            //break;
            a_j = mean_step;
        }
        
        if (isnan(a_j) || isinf(a_j))
        {
            fprintf(stderr,"nan/inf value\n");
            success = 0;
            break;
        }

        /* evaluate function at a_lo (don't need dir. derivative) */
        for (j=0;j<N;j++)
        {
            x_update[j] = x[j] + a_lo * p[j];
        } 
        phi_a_lo = misfit(x_update, N, aux);


        /* evaluate function at a_j. directional derivative will be         */
        /* calculated later if needed                                       */
        for (j=0;j<N;j++)
        {
            x_update[j] = x[j] + a_j * p[j];
        } 
        phi_aj = misfit(x_update, N, aux);

        printf("Trial step length %E misfit: %E\n", a_j, phi_aj);

        if (phi_aj > phi0 + WOLFE_C1 * a_j * phi_prime0
            || phi_aj >= phi_a_lo)
        {
            /* trial step length violates Armijo or                         */
            /* trial step evaluate to a function value greater than a_lo    */
            /* need to reduce step length                                   */
            a_hi = a_j;
        }
        else 
        {
            /* trial step length satisfies Armijo and function value is     */
            /* less than a_lo                                               */
            /* check curvature                                              */
            grad(x_update, gradient, N, aux);
            phi_prime_aj = vec_dot_prod(gradient, p, N);
            printf("Trial step length %E dir deriv: %E\n", a_j, phi_prime_aj);
            if (fabs(phi_prime_aj) <= -WOLFE_C2 * phi_prime0)
            {
                /* trial step length satisfies strong curvature condition   */
                /* accept trial step length                                 */
                bb_star = a_j;
                success = 1;
                break;
            }

            if (phi_prime_aj * (a_hi - a_lo) >= 0)
            {
                a_hi = a_lo;
            }

            /* satisfies Armijo, but not curvature                          */
            /* increase step length                                         */
            a_lo = a_j;


        }
    }

    free(x_update);
    free(gradient);

    if (!success)
    {
        fprintf(stderr,"Line search zoom failed to return an "
        "appropriate step length\n");
        //exit(-90);
        return -1;
    }

    return bb_star;
}


/*****************************************************************************/

double ls_cubic_interpolation(const double *x, 
               const double *p,
               const int N,
               optimization_parameters *opt_par,
               const double bound1,
               const double bound2,
               void *aux
               )
{
    /* calculates a new trial step length from cubic interpolation  */
    /* interpolation is based on two function values and the        */
    /* corresponding directional derivatives                        */
    /*                                                              */
    /* This interpolation requires two function evaluations and     */
    /* two gradient evaluations. If the function or gradient is     */
    /* expensive to evaluate, there are probably better ways to     */
    /* calculate an appropriate step length                         */
    /*                                                              */
    /* see Nocedal and Wright 2006, equation 3.59                   */
    /*                                                              */
    /* inputs:                                                      */ 
    /*      x: current point in parameter space                     */
    /*      p: search direction                                     */
    /*      N: number of dimensions                                 */
    /*      misfit: function that takes point and returns objective */
    /*              function value. x and N as defined above        */
    /*      grad: function that takes point and returns gradient    */
    /*            value (returned as parameter)                     */
    /*      bound1, bound2: bounds of step length search            */
    /*                                                              */
    /*  new trial step length is returned                           */
    
    
    double (*misfit)(const double *, const int, void *) = opt_par->objective;
    void (*grad)(const double *, double *, int, void *) = opt_par->grad;


    double *gradient = opt_allocate(sizeof(*gradient) * N);
    double *x_update = opt_allocate(sizeof(*x_update) * N);

    double f1, f2, g1, g2;

    int i;

    /* evaluate objective function and directional derivative at bound1 */
    for (i=0;i<N;i++)
    {
        x_update[i] = x[i] + bound1 * p[i];
    } 
    f1 = misfit(x_update, N, aux);
    grad(x_update, gradient, N, aux);
    g1 = vec_dot_prod(gradient, p, N);

    printf("misfit at    %E: %E\n", bound1, f1);
    printf("dir deriv at %E: %E\n", bound1, g1);

    /* evaluate objective function and directional derivative at bound2 */
    for (i=0;i<N;i++)
    {
        x_update[i] = x[i] + bound2 * p[i];
    } 
    f2 = misfit(x_update, N, aux);
    grad(x_update, gradient, N, aux);
    g2 = vec_dot_prod(gradient, p, N);


    printf("misfit at    %E: %E\n", bound2, f2);
    printf("dir deriv at %E: %E\n", bound2, g2);

    /*calculate cubic approximateion */
    double d1 = g1 + g2 - 3 * ((f1 - f2)/(bound1 - bound2));

    /* note that both bound1 and bound2 should always be positive   */
    /* they should never be equal                                   */
    /* d2_sign is sign(bound2 - bound1)                             */
    double d2_sign = bound2 > bound1 ? 1.0 : -1.0;
    if ((bound2 - bound1) > 0)
        d2_sign = 1.0;
    else if ((bound2 - bound1) < 0)
        d2_sign = -1.0;
    else
        d2_sign = 0.0;


    double d2 = d2_sign * sqrt(d1 * d1 - g1 * g2);

    double a_new = bound2 - (bound2 - bound1) * ((g2 + d2 - d1)/(g2 - g1 + 2 * d2));

    /* the minimum is always either at a_new as calculated above orj    */
    /* at one of the end points. Check those here.                      */
    //a_new = (a_new < bound1) ? a_new : bound1;
    //a_new = (a_new < bound2) ? a_new : bound2;

    free(x_update);
    free(gradient);
    return a_new;
}


/*****************************************************************************/


double ls_bracket(const double *x, 
               const double *p,
               const int N,
               optimization_parameters *opt_par,
               double alpha0,
               void *aux
               )
{

    double (*misfit)(const double *, const int, void *) = opt_par->objective;
    void (*grad)(const double *, double *, int, void *) = opt_par->grad;


    int success = 0;
    double b_star;

    int steplength_max = 50;
    double beta_max = 5e30;
    double beta0 = 0;
    double beta1 = alpha0; /* initial step length */

    double *gradient = opt_allocate(sizeof(*gradient) * N);
    double *x_update = opt_allocate(sizeof(*x_update) * N);


    double p0;      /* phi(0)  */
    double pp0;     /* phi'(0) */
    /* evaluate objective and dir. derivative at initial point */
    p0 = misfit(x, N, aux);
    grad(x, gradient, N, aux);
    pp0 = vec_dot_prod(gradient, p, N);


    double pi;              /* phi(i)       */
    double ppi;             /* phi'(i)      */
    double pi_1  = p0;      /* phi(i-1)     */
    //double ppi_1 = pp0;     /* phi'(i-1)    */
    double b_i   = beta1;   /* beta_i       */
    double b_i_1 = beta0;   /* beta_{i-1}   */
    

    int i, j;
    for (i=0;i<steplength_max;i++)
    {
        printf("Starting bracketing, iteration %d\n",i);
        printf("Bracket bounds: %E, %E\n",b_i_1, b_i);
        fflush(stdout);

        /* calculate objective and derivative values at b_i */
        for (j=0;j<N;j++)
        {
            x_update[j] = x[j] + b_i * p[j];
        }
        pi = misfit(x_update, N, aux);

        if ((pi > p0 + WOLFE_C1 * b_i * pp0) || (pi >= pi_1 && i>0))
        {
            b_star = ls_zoom(x, p, N, opt_par, b_i_1, b_i, aux);
            if (b_star < 0)
            {
                /* ls_zoom failed - try bigger steps */
                b_i = 2 * alpha0;
                b_i_1 = alpha0;
                continue;
            }
            success = 1;
            break;
        }

        grad(x_update, gradient, N, aux);
        ppi = vec_dot_prod(gradient, p, N);
        if (fabs(ppi) <= -WOLFE_C2 * pp0)
        {
            b_star = b_i;
            success = 1;
            break;
        }

        if (ppi >= 0.0)
        {
            b_star = ls_zoom(x, p, N, opt_par, b_i, b_i_1, aux);
            if (b_star < 0)
            {
                /* ls_zoom failed - try bigger steps */
                b_i = 2 * alpha0;
                b_i_1 = alpha0;
                continue;
            }
            success = 1;
            break;
        }

        /* update b_i and b_i_1 */
        b_i_1 = b_i;
        b_i = 2 * b_i;

        if (b_i > beta_max)
        {
            fprintf(stderr,"Error: ls_bracket: max step length reached\n");
            success = 0;
            break;
        }
    }


    free(x_update);
    free(gradient);


    if (!success)
    {
        fprintf(stderr,"Line search bracket failed to return an "
        "appropriate step length\n");
        exit(-90);
    }

    printf("Step length successfully selected. Step length: %E\n",b_star);
    fflush(stdout);
    return b_star;
}




/*****************************************************************************/

int check_angle_restart(double *grad, double *p, int N)
{
    int restart = 0;
    double angle_tol_deg = 95.0;
    double angle_tol_rad = angle_tol_deg * (M_PI / 180.0); 

    double grad_norm = vec_norm(grad, N);
    double p_norm = vec_norm(p, N);
    double pg = vec_dot_prod(p, grad, N);
    double ang = pg / (p_norm * grad_norm);
    if(ang > cos(angle_tol_rad))
    {
        restart = 1;
    }

    return restart;
}


/*****************************************************************************/


double back_track_quad(double phi0, double phip0, double phi_a0, double a0)
{
    /* see Nodedal and Wright, p. 58 (equation 3.57 and 3.58) */

    double a1;

    double num = phip0 * a0 * a0;
    double den = phi_a0 - phi0 - phip0 * a0;

    a1 = -0.5 * (num / den);
    return a1;
}


double back_track_cubic(double phi0, 
                        double phip0,
                        double phi_a0,
                        double phi_a1,
                        double a0,
                        double a1)
{
    /* see Nocedal and Wright 2006, p. 58 */


    double a, b, a2;

    double a02 = a0 * a0;   /* a0 squared   */
    double a12 = a1 * a1;   /* a1 squared   */
    double a03 = a02 * a0;  /* a0 cubed     */
    double a13 = a12 * a1;  /* a1 cubed     */ 


    double v1 = phi_a1 - phi0 - phip0 * a1;
    double v2 = phi_a0 - phi0 - phip0 * a0;

    double norm = 1.0 / (a02 * a12 * (a1 - a0));


    a = norm * (a02 * v1 - a12 * v2);
    b = norm * (-a03 * v1 + a13 * v2);


    a2 = (-b + sqrt(b * b - 3 * a * phip0)) / (3 * a);

    return a2;
}




/*****************************************************************************/


double back_track(double *x,
                  double *p,
                  int N,
                  optimization_parameters *opt_par,
                  double alpha,
                  void *aux
                  )
{
    static double a[2]; /* holds a set of step lengths that can be used for bracketing */
    a[0] = 0;
    a[1] = 0;


    double (*misfit)(const double *, const int, void *) = opt_par->objective;
    void (*grad)(const double *, double *, int, void *) = opt_par->grad;

    int i;

    double phi0, phip0;
    double a0, a1, a2;
    double phi_a0, phi_a1, phi_a2;

    double *gradient = opt_allocate(sizeof(*gradient) * N);
    double *x_update = opt_allocate(sizeof(*x_update) * N);

    /* evaluate initial misfit and directional derivative */
    phi0 = misfit(x, N, aux);
    grad(x, gradient, N, aux);
    phip0 = vec_dot_prod(gradient, p, N);


    /* try a0*/
    a0 = alpha;
    a[1] = a0;
    for (i=0;i<N;i++)
    {
        x_update[i] = x[i] + a0 * p[i];
    }
    phi_a0 = misfit(x_update, N, aux);

    /* check Armijo condition */
    if (phi_a0 <= phi0 + WOLFE_C1 * a0 * phip0)
    {
        /* condition satisfied */
        free(x_update);
        free(gradient);
        return a0;
    }

    /* try quadratic approximation */
    a1 = back_track_quad(phi0, phip0, phi_a0, a0);
    a[0] = a1;

    for (i=0;i<N;i++)
    {
        x_update[i] = x[i] + a1 * p[i];
    }
    phi_a1 = misfit(x_update, N, aux);

    /* check Armijo condition */
    if (phi_a1 <= phi0 + WOLFE_C1 * a1 * phip0)
    {
        /* condition satisfied */
        free(x_update);
        free(gradient);
        return a1;
    }

    /* move into cubic approximations */
    int cubic_max = 10; /* <----------May need to adjust max number of iterations */
    int j;

    for (j=0;j<cubic_max;j++)
    {
        a2 = back_track_cubic(phi0, phip0, phi_a0, phi_a1, a0, a1);

        /* make sure a2 is not too close to a1 */
        if (fabs(a2 - a1) / a1 < 0.01) /* <--------- May need to adjust the factor of 0.01 here */
            a2 = 0.5 * a1;

        /* make sure a2 is not too much smaller than a1 */
        if (a2 < (0.01 * a1)) /* <-------- May need to adjust the factor of 0.01 here */
            a2 = 0.5 * a1;

        a[1] = a[0];
        a[0] = a2;


        for (i=0;i<N;i++)
        {
            x_update[i] = x[i] + a2 * p[i];
        }
        phi_a2 = misfit(x_update, N, aux);

        /* check Armijo condition */
        if (phi_a2 <= phi0 + WOLFE_C1 * a2 * phip0)
        {
            /* condition satisfied */
            break;
        }

        /* keep only the two newest values */
        a0 = a1;
        phi_a0 = phi_a1;
        a1 = a2;
        phi_a1 = phi_a2;
    }

    
    if (j >= cubic_max)
    {
        fprintf(stderr,"Backtrack failed to return an "
                       "appropriate step length\n");
        //exit(1);
    }



    free(x_update);
    free(gradient);
    return a2;
}










/*****************************************************************************/

double max3(double x1, double x2, double x3)
{
    double mm = (x1 > x2) ? x1 : x2;
    mm = (mm > x3) ? mm : x3;
    return mm;

}



void modified_chol(double *A, int N)
{

    /* make temp copy of A */
    double *temp = opt_allocate(sizeof(*temp) * N * N);
    double *c = opt_allocate(sizeof(*c) * N * N);
    double *d = opt_allocate(sizeof(*d) * N);
    memcpy(temp, A, sizeof(*temp) * N *N);


    int i, j, s;

    double sum = 0;

    double beta = 1.0;
    double delta = 1.0e-3;
    double theta = 0.0;

    for (j = 0; j < N; j++)
    {
        for (s=0;s<j;s++)
        {
            sum = d[s] * A[j+s*N] * A[j+s*N];
        }

        c[j+j*N] = temp[j+j*N] - sum; 
        A[j+j*N] = 1.0;

        //d[j] = c[j+j*N];
        d[j] = max3(fabs(c[j+j*N]),pow(theta / beta,2.0),delta);
        theta = DBL_MIN;
        for (i = j+1; i<N;i++)
        {
            for (s=0;s<j;s++)
            {
                sum = d[s] * A[i+s*N] * A[j+s*N];
            }

            c[i+j*N] = temp[j+i*N] -sum;
            theta = theta > c[i+j*N] ? theta : c[i+j*N];
            A[i+j*N] = c[i+j*N] / d[j];
        }
    }

    for (j=0;j<N;j++)
    {
        for (i=j;i<N;i++)
        {
            A[i+j*N] = A[i+j*N] * sqrt(d[j]);
        }
    }

    //printf("Diag:\n");
    //for (i=0;i<N;i++)
    //    printf("%f\n",d[i]);


    free(d);
    free(c);
    free(temp);
}




void* opt_allocate(long nbytes)
{
    void *ptr = malloc(nbytes);
    if (!ptr)
    {
        fprintf(stderr, "Error: optimization.c : memory allocation failure\n");
        exit (1);
    }
    return ptr;
}