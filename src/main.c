#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "fwi.h"

#include "ibool_model.h"
#include "mask_kernels.h"
#include "unique_model.h"

#include "bilinear_interp.h"
#include "gauss_conv.h"

#include "optimization.h"

#include "read_parameters.h"

#include "fwi_misfit.h"

int main(int argc, char **argv)
{
   if (argc != 2)
   {
      fprintf(stderr,"Please specify a parameter file\n");
      exit(1);
   }


   printf("Starting FWI program\n");

   fwi ff_;
   fwi *ff = &ff_;

   parameters par_;
   parameters *par = &par_;

   /* read parameters */
   printf("Reading parameters...");
   read_parameters(par,argv[1]);
   printf("done\n");

   /* print parameters to stdout */
   print_parameters(par);


   /* initialize all models and data - 
   this also reads in initial model(ibool) */ 
   printf("Initializing fwi...");
   initialize_fwi(ff,par);
   printf("done\n");
   
   /* convert ibool model to unique */
   printf("Converting model...");
   ibool2unique(ff->imod, ff->umod, ff->i2u);
   printf("done\n");
   model_min_max(ff->umod);
   
   /* convert model from vp/vs to rmk */
   printf("Converting to rmk and nat log...");
   vpvs2rmk(ff->umod);
   ln_mod(ff->umod);
   printf("done\n");

   /* interpolate unique model to regular grid - hold interpolation
      structure in memory */
   printf("Interpolating model...");
   interpolate_model(ff->umod, ff->umod_reg, ff->sim2reg);
   printf("done\n");
   printf("Interpolated Model: ");
   model_min_max(ff->umod_reg);

   /* spatial filter interpolated model */
   printf("Applying spatial filter...");
   gauss_conv_model(ff->umod_reg, ff->umod_filt, par->sp_lp_cutoff);
   printf("done\n");

   /* downsample interpolated model */
   printf("Downsampling model...");
   int dec_factor = (int)(par->dx_inv / par->dx_reg);
   decimate_model(ff->umod_inv,ff->umod_filt,dec_factor);
   printf("done\n");
   fflush(stdout);


   /* DONE WITH PRELIM */
   /*************************************************************************/

   /* assign optimization parameters */

   optimization_parameters opt_par_;
   optimization_parameters *opt_par = &opt_par_;

   opt_par->max_iterations = par->max_iterations;
   opt_par->stopping_tolerance = par->stopping_tol;
   opt_par->initial_step_length = par->init_step_mult;
   
   opt_par->sd = par->search_direction;
   opt_par->LBFGS_mem = par->LBFGS_mem;
   opt_par->angle_restart_allowed = par->angle_restart;
   
   opt_par->backtrack_only = par->backtrack_only;
   opt_par->bracket_only = par->bracket_only;
   opt_par->bracket_zoom_iterations = 20;
   opt_par->bracket_zoom_tolerance = 1e-14;

   opt_par->objective = fwi_misfit;
   opt_par->grad = fwi_gradient;
   opt_par->hess = NULL;

   opt_par->strong_curvature = par->strong_curvature;
   opt_par->wolfe_c1 = par->wolfe_c1;
   opt_par->wolfe_c2 = par->wolfe_c2;


   /* perform optimization */

   /* convert initial model to a vector */
   int N = get_ndim(ff);
   double *x_in   = malloc(sizeof(*x_in) * N);
   double *x_out  = malloc(sizeof(*x_out) * N); 
   if (!x_in || !x_out)
   {
      fprintf(stderr,"Error: main: memory allocation failure\n");
      exit(1);
   }

   model2vector(ff, x_in);

   write_new_model(ff,x_in);

   line_search(x_in,x_out, N, opt_par,ff);


   /* write the final model to a file */
   FILE *mod_out = fopen("./model_out.csv","w");

   exp_mod(ff->umod_inv);
   rmk2vpvs(ff->umod_inv);


   int numele = ff->umod_inv->num_ele;
   int i;
   double xx, zz, vp, vs, rho;
   fprintf(mod_out, "X, Z, DENS, VP, VS\n");
   for (i=0;i<numele;i++)
   {
      xx = ff->umod_inv->x[i];
      zz = ff->umod_inv->z[i];
      vp = ff->umod_inv->vp[i];
      vs = ff->umod_inv->vs[i];
      rho = ff->umod_inv->dens[i];

      fprintf(mod_out, "%E, %E, %E, %E, %E\n", xx, zz, rho, vp, vs);
   }
   fclose(mod_out);

   free(x_in);
   free(x_out);



   free_fwi(ff);


   return 0;
}
