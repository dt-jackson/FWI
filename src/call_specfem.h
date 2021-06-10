#ifndef CALL_SPECFEM_H
#define CALL_SPECFEM_H

/* success exit code */
#define SF_EXIT_SUCCESS 0

/* wrapper functions to call specfem */
void specfem2d_forward(const char *sf_forward_bin);
void specfem2d_adjoint(const char *sf_adjoint_bin);


#endif