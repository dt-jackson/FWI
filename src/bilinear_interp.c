#include <float.h>
#include <libqhull_r/geom_r.h>
#include <libqhull_r/libqhull_r.h>
#include <libqhull_r/mem_r.h>
#include <libqhull_r/poly_r.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <time.h>

#include <libqhull_r/qhull_ra.h>

#include "bilinear_interp.h"


#ifndef LOOP_MAX
#define LOOP_MAX 500
#endif

#ifndef TOLERANCE
#define TOLERANCE 0.00005
#endif

/* The functions here perform bilinear interpolation using a Delaunay        */
/* triangulation of the points                                               */
/* An interpolation structure must first be built from the triangulation,    */
/* but once that is calculated, subsequent interpolations can be performed   */
/* quickly                                                                   */
/*                                                                           */
/* This must be linked against the reentrant version of libqhull             */
/* i.e., libqhull_r                                                          */
/*                                                                           */
/*****************************************************************************/

void bilinear_interp(double  *restrict r, double *restrict r_i, int nr_i, interp2d_t *restrict pp)
{
    /* interpolates from r to r_i                       */
    /* nr_i is the number of points in r_i              */
    /* pp holds the information about the interpolation */
    /* call bilinear_setup to generate pp               */
    
    int i;
    for (i=0;i<nr_i;i++)
    {
        r_i[i] = 
          r[pp[i].vert_list[0]] * pp[i].lam[0]
        + r[pp[i].vert_list[1]] * pp[i].lam[1]
        + r[pp[i].vert_list[2]] * pp[i].lam[2];
    }

} 



/*****************************************************************************/


void bilinear_setup(double *x,double *y,double *qx,double *qy,int n,int nq, interp2d_t *pp)
{
    /* bilinear interpolation (ndims = 2) */
    /* note that this function just sets up the interpolation structure,    */
    /* but doesn't actually perform the interpolation                       */
    int ndim = 2;

    int i;

    /* copy data into the format qhull wants                  */
    /* each qpoint gets 1 extra dimension so it can be lifted */
    /* because of the way we're calling qhull, we don't need  */
    /* to allocate an extra dimension to points               */
    coordT *points = malloc((ndim)*n*sizeof(*points));
    coordT *qpoints = malloc((ndim+1)*nq*sizeof(*qpoints));
    if (points == NULL || qpoints == NULL)
    {
        fprintf(stderr,"Memory allocation failure in bilinear_setup\n");
        exit(1);
    }

    for (i=0;i<n;i++)
    {
        points[i*(ndim)] = x[i];
        points[i*(ndim)+1] = y[i];

    }

    for (i=0;i<nq;i++)
    {
        qpoints[i*(ndim+1)] = pp[i].pt[0] = qx[i];
        qpoints[i*(ndim+1)+1] = pp[i].pt[1] = qy[i];
        /* this line "lifts" the points */
        qpoints[i*(ndim+1)+2] = qx[i] * qx[i] + qy[i] * qy[i];
    }

    /* declare qhull variables */
    int exitcode;        
    boolT ismalloc = False;
    char flags[250];
    sprintf(flags,"qhull d QJ");

    qhT qh_qh;
    qhT *qh = &qh_qh;
    qh_zero(qh,stderr);
    facetT *ff;
    vertexT *vertex, **vertexp;

    /* call qhull to get delaunay triangulation */    
    #ifdef FWIVERBOSE
    printf("Calculating Delaunay Triangulation...\n");
    #endif

    exitcode = qh_new_qhull(qh,ndim,n,points,ismalloc,flags,NULL,stderr);
    
    #ifdef FWIVERBOSE
    printf("Triangulation complete\n");
    printf("Locating points...\n");
    #endif




    /* find out which triangle each query point is in */
    facetT *facet = NULL;
    for (i=0;i<nq;i++)
    {
        //ff = findDelaunay2d(qh,  qpoints+i*(ndim+1),pp+i);
        //if (!ff)
        //    fprintf(stderr,"Could not find triangle\n");
        facet = find_tri(qh,  qpoints+i*(ndim+1),pp+i,facet);
    }

    #ifdef FWIVERBOSE
    printf("All points successfully located\n");
    #endif

    /* free qhull structures */
    /* memory remaining after qh_memfreeshort, used if !qh_NOmem  */
    int curlong, totlong;     
    #ifdef qh_NOmem
        qh_freeqhull(qh, qh_ALL);
    #else
        /* free long memory  */
        qh_freeqhull(qh, !qh_ALL);
        /* free short memory and memory allocator */
        qh_memfreeshort(qh, &curlong, &totlong);
        if (curlong || totlong)
            fprintf(stderr, "qhull internal warning (bilinear_setup): "
                "did not free %d bytes of long memory (%d pieces)\n", 
                totlong, curlong);
    #endif

    /* free other arrays     */
    free(qpoints);
    free(points);
}


/*****************************************************************************/



/*****************************************************************************/


/* helper functions to find the max/min index of a 3 element array */
int max_elem(double* x)
{
   int k = x[0] > x[1] ? 0 : 1; 
   return x[k] > x[2] ? k : 2;
}

int min_elem(double *x, int k)
{
    double m = x[0];
    int mi = 0;
    int i;
    for (i=0;i<k;i++)
    {
        if (x[i] < m)
        {
            m = x[i];
            mi = i;
        }
    }
    return mi;

    //int k = x[0] < x[1] ? 0  : 1;
    //return x[k] < x[2] ? k : 2;
}

/*****************************************************************************/
/* function to determine which Delaynay triangle a point is in               */
/* tries a fast way first, but then moves to an exhaustive search if no      */
/* triangle can be found                                                     */
/* if no triangle can be found after the exhaustive search, the points are   */
/* joggled and a new exhaustive search is started. If this still fails, the  */
/* program exits with an error                                               */
/* see Mucke et al. 1996. Fast randomized point location without             */
/* preprocessing in two- and three-dimensional Delaunay triangulations       */ 

facetT* find_tri(qhT *qh,coordT *point, interp2d_t *pp, facetT *facet)
{
    int dim = 2;
    boolT isoutside;
    realT bestdist;

    /* use qhull to find the first facet */ 
    if (facet == NULL)
    {
        /* note this only happens the first time or if facet is */
        /* reset to NULL                                        */
        facet= qh_findbestfacet(qh, point, qh_ALL, &bestdist, &isoutside);
    }

    /* check if we are in the triangle from the previous */
    /* function call                                     */
    int good = in_triangle(qh,facet,point,pp);

    /* if we found it already, return */
    if (good && !facet->upperdelaunay)
        return facet;

    good = 0; 
    /* search neighboring facets */
    double *fcenter,*ncenter;
    double fvec[3], nvec[2];
    double dd[10];
    facetT *neighbor, **neighborp;
    vertexT *vertex, **vertexp;
    int vert_count;

    fcenter = qh_getcentrum(qh, facet);
    double p0[2];
    p0[0] = fcenter[0];
    p0[1] = fcenter[1];
    //printf("starting point: %f, %f\n",p0[0],p0[1]);

    double p1[2];
    p1[0] = point[0];
    p1[1] = point[1];
    //printf("target point: %f, %f\n",p1[0],p1[1]);

    double xx = p1[0] - p0[0];
    double yy = p1[1] - p0[1];

    double xxx;
    double yyy;
    int sign[3];
    int pass_through;

    #ifdef EXHAUSTIVE_INTERP
    int loopcount = LOOP_MAX;
    #else
    int loopcount = 0;
    #endif
    facetT *old_facet = facet;
    double normal_tol = 0.0001;
    while (!good && facet && loopcount < LOOP_MAX)
    {
        //printf("facet: ");
        //FOREACHvertex_(facet->vertices)
        //{
        //    printf("%f, %f\t",vertex->point[0],vertex->point[1]);
        //}
        //printf("\n");
        //printf("is upper: %u\n",facet->upperdelaunay);
        //printf("%f, %f, %f", facet->normal[0], facet->normal[1], facet->normal[2]);
        //printf("\n");

        int k=0;
        FOREACHneighbor_(facet)
        {
            if (!neighbor->upperdelaunay)
            {   
                int kk=0; 
                FOREACHvertex_(neighbor->vertices)
                {
                    xxx = vertex->point[0];
                    yyy = vertex->point[1];

                    double val = xx*(yyy-p0[1]) - yy*(xxx-p0[0]);
                    if (val > 0) sign[kk] = 1;
                    if (val < 0) sign[kk] = -1;
                    kk++;
                }

                fcenter = qh_getcentrum(qh, neighbor);
                dd[k]  = (point[0] - fcenter[0]) * (point[0] - fcenter[0]);
                dd[k] += (point[1] - fcenter[1]) * (point[1] - fcenter[1]);
                
                if ((sign[0] == sign[1] && sign[0] == sign[2]) || neighbor == old_facet)
                    dd[k] = 1e10;

                //fcenter = qh_getcentrum(qh, neighbor);
                //dd[k]  = (point[0] - fcenter[0]) * (point[0] - fcenter[0]);
                //dd[k] += (point[1] - fcenter[1]) * (point[1] - fcenter[1]);

                //dd[k] += (point[2] - fcenter[2]) * (point[2] - fcenter[2]);
                //qh_distplane(qh,point,neighbor,dd+k); 
                ///* distances >0 are above */
                //if (dd[k]>0) dd[k]=0;
                //dd[k] *= -1;
                k++;
            }
        }
        if (k==0)
        {
            loopcount = LOOP_MAX;
            //#ifdef VERBOSE
            //printf("no neighbor\n");
            //#endif
            break;
        }

        /* find min distance and go there */
        k = min_elem(dd,k);
        old_facet = facet;
        facet = (facetT*)facet->neighbors->e[k].p;
        good = in_triangle(qh,facet,point,pp);

        /* for some reason, some upper dealunay facets still get to this point */
        /* eliminate them here                                                 */
        if (facet->upperdelaunay == 1)
            good = 0;

        loopcount++;
    }

    /* if we couldn't find the triangle before, start an */
    /* exhaustive search                                 */
    if (loopcount == LOOP_MAX && !good)
    {
        #ifdef VERBOSE
        printf("starting exhaustive search\n");
        printf("%f, %f\n",point[0],point[1]);
        #endif

        FORALLfacets{
            if (!facet->upperdelaunay)
            {
                if (in_triangle(qh, facet, point, pp))
                {
                    return facet;
                }
            }
        }
        #ifdef VERBOSE
        printf("facet not found for point ");
        printf("%f, %f\n",point[0],point[1]);
        #endif

        /* if we haven't found it yet, set facet to null */
        /* to try next thing                             */
        facet = NULL;
    }

    /* if triangle is still not found, it is outside of the */
    /* convex hull, possibly due to floating point error    */  
    /* joggle the points a little and try again             */
    if (!facet)
    {
        #ifdef VERBOSE
        printf("Joggling points\n");
        #endif
        int llcc=0;
        double mult = 0.5 * FLT_EPSILON;
        srand((unsigned int)time(NULL));
        while (!facet)
        {
            /* amount points are joggled is doubled each iteration */
            mult *=2;

            /* random values to add or subtract */
            double j1 = mult * (((double)rand() / (double)RAND_MAX) - 0.5);
            double j2 = mult * (((double)rand() / (double)RAND_MAX) - 0.5);


            /* since we don't know which way we need to go, */
            /* try every combination of adding/subtracting  */
            double tmp_pt[12];

            tmp_pt[0] = point[0] + j1;
            tmp_pt[1] = point[1] + j2;
            tmp_pt[2] = tmp_pt[0] * tmp_pt[0] + tmp_pt[1] * tmp_pt[1];

            tmp_pt[3] = point[0] - j1;
            tmp_pt[4] = point[1] + j2;
            tmp_pt[5] = tmp_pt[3] * tmp_pt[3] + tmp_pt[4] * tmp_pt[4];

            tmp_pt[6] = point[0] + j1;
            tmp_pt[7] = point[1] - j2;
            tmp_pt[8] = tmp_pt[6] * tmp_pt[6] + tmp_pt[7] * tmp_pt[7];

            tmp_pt[9] = point[0] - j1;
            tmp_pt[10] = point[1] - j2;
            tmp_pt[11] = tmp_pt[9] * tmp_pt[9] + tmp_pt[10] * tmp_pt[10];

            FORALLfacets{
                if (!facet->upperdelaunay)
                {
                    for (int i=0;i<4;i++)
                    {
                        if (in_triangle(qh,facet,tmp_pt+i*3,pp))
                        {
                            #ifdef VERBOSE
                            printf("facet found for point ");
                            printf("%f, %f\n",tmp_pt[0+i*3],tmp_pt[1+i*3]);
                            printf("It took %d iterations\n",llcc);
                            printf("Found in direction %d\n",i);
                            //printf("%.25f\n",FLT_EPSILON);
                            #endif

                            return facet;
                        }
                    }
                }
            }
            facet = NULL;
            llcc++;
            if (mult > TOLERANCE)
                break;
        }
    }
    if (facet == NULL)
    {
        fprintf(stderr,"Unable to find Delaunay region for point ");
        fprintf(stderr,"%f, %f\n",point[0],point[1]);
        fprintf(stderr,"Exiting with error code 60\n");
        exit(60);
    }

    if(facet->upperdelaunay)
    {
        fprintf(stderr,"upper delaunay\n");
        fprintf(stderr, "%f, %f\n", point[0], point[1]);
    }
    
    return facet;
}

/*****************************************************************************/
/* calculate barycentric weights and determine if a point is in a given      */
/* triangle                                                                  */


int in_triangle(qhT *qh, facetT *ff,coordT *qp,interp2d_t *pp)
{
    /* check if qp is in facet ff and calculate barycentric weights */
    /* returns 1 if qp is in triangle and 0 otherwise               */

    /* quick return if ff is null */
    if (ff == NULL)
        return 0;
    
    int k = 0;

    /* loop through each vertex and get info */
    vertexT *vertex, **vertexp;
    FOREACHvertex_(ff->vertices)
    {
        pp->x[k] = vertex->point[0];
        pp->y[k] = vertex->point[1];
        pp->vert_list[k] = qh_pointid(qh,vertex->point);
        k++;
    }

    /* calculate barycentric coordinates */
    barycentric_calc2d(pp->x, pp->y, qp, pp->lam);

    /* if any barycentric coords are negative qp is outside */
    /* of triangle - return 0                               */
    if (pp->lam[0] < -TOLERANCE || pp->lam[1] < -TOLERANCE || pp->lam[2] < -TOLERANCE)
        return 0;
    //if (pp->lam[0] < 0 || pp->lam[1] < 0 || pp->lam[2] < 0)
    //    return 0;

    /* return 1 if all weights were >=0 */
    return 1;
}

/*****************************************************************************/
/* calculate barycentric weights for a given triangle and point              */

void barycentric_calc2d(double *x,double *y,double *qp,double *lam)
{
    /* calculate barycentric weights for qp */
    /* result returned in lam               */

    /* inverse of denominator */
    double den = 1 / ((y[1] - y[2]) * (x[0] - x[2]) 
               + (x[2] - x[1]) * (y[0] - y[2]));
    
    lam[0] = (y[1]-y[2]) * (qp[0]-x[2]) + (x[2] - x[1]) * (qp[1] - y[2]);
    lam[1] = (y[2]-y[0]) * (qp[0]-x[2]) + (x[0] - x[2]) * (qp[1] - y[2]);
    lam[0] *= den;
    lam[1] *= den;
    lam[2] = 1.0 - lam[0] - lam[1];
}

/*****************************************************************************/


void write_interp_struct2d(const interp2d_t *interp, const int N, const char *fname)
{
    FILE *fid = fopen(fname, "wb");
    if (!fid)
    {
        fprintf(stderr,"Failure opening %s for writing\n",fname);
        exit(1);
    }

    int i;
    for(i=0;i<N;i++)
    {
        /*  write point  list */
        fwrite(interp[i].vert_list,sizeof(*interp[i].vert_list), 3, fid);

        /* write weights */
        fwrite(interp[i].lam,sizeof(*interp[i].lam), 3, fid);
    }

    fclose(fid);
}

void read_interp_struct2d(interp2d_t *interp, const int N, const char *fname)
{
    FILE *fid = fopen(fname, "rb");
    if (!fid)
    {
        fprintf(stderr,"Failure opening %s for reading\n",fname);
        exit(1);
    }

    int i;
    for(i=0;i<N;i++)
    {
        /*  read point  list */
        fread(interp[i].vert_list, sizeof(*interp[i].vert_list), 3, fid);

        /* read weights */
        fread(interp[i].lam,sizeof(*interp[i].lam), 3, fid);
    }

    fclose(fid);
}