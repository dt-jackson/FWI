#include <stdlib.h>
#include <stdio.h>

#include "call_specfem.h"


void specfem2d_forward(const char *sf_forward_bin)
{
    int status = system(sf_forward_bin);
    if (status != SF_EXIT_SUCCESS)
    {
        fprintf(stderr, "SPECFEM failure: %s\n",sf_forward_bin);
        exit(status);
    }
}


void specfem2d_adjoint(const char *sf_adjoint_bin)
{
    int status = system(sf_adjoint_bin);
    if (status != SF_EXIT_SUCCESS)
    {
        fprintf(stderr, "SPECFEM failure: %s\n",sf_adjoint_bin);
        exit(status);
    }
}