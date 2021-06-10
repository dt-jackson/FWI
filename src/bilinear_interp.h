#ifndef BILINEAR_INTERP_H
#define BILINEAR_INTERP_H

#include <libqhull_r/qhull_ra.h>


struct interp2d
{
    int nq;             /* number of query points          */
    double pt[2];       /* x,y coordinate of query point   */
    int vert_list[3];   /* verticies of triangle around pt */
    double x[3];        /* x coords of verticies           */
    double y[3];        /* y coords of verticies           */
    double lam[3];      /* barycentric weights             */
};
typedef struct interp2d interp2d_t;

void bilinear_setup(double *x,double *y,double *qx,double *qy,int n,int nq, interp2d_t *pp);
void bilinear_interp(double  *r, double *r_i, int nr_i, interp2d_t *pp);

facetT* find_tri(qhT *qh,coordT *point, interp2d_t *pp,facetT *facet);
int in_triangle(qhT *qh, facetT *ff,coordT *qp, interp2d_t *pp);
void barycentric_calc2d(double *x,double *y,double *qp,double *lam);

//int max_elem(double* x);
//int min_elem(double *x);

#endif