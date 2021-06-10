#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "stats.h"

#include "median_filt.h"

void median_filter1d(double *x, int N, int filt_size)
{
    if (filt_size <=1) return;

    /* make a temp copy of x */
    double *y = malloc(sizeof(*y) * N);
    if (!y)
    {
        fprintf(stderr, "Error: median_filter1d: memory allocation failure\n");
        exit(1);
    }


    memcpy(y, x, sizeof(*y) * N);


    int half_size = filt_size/2;


    double (*med)(const double *, const int);
    
    if (filt_size <= 9)
        med = med_small;
    else
    {
        //med = median_of_medians;
        med = median;
    }


    int i;

    double *p = y;
    /* handle first filt_size/2 elements */
    for (i=0;i<half_size;i++)
    {
        x[i] = med(p, i + 1 + half_size);
    }
    
    /*loop through each element in the middle */
    for (i=half_size;i< N - half_size;i++, p++)
    {
        x[i] = med(p, filt_size);
    }

    /* handle last filt_size/2 elements */
    int c = 1;
    for (i= N-half_size; i<N; i++, p++)
    {
        x[i] = med(p, filt_size -c);
        c++;
    }

    free(y);
}


void median_filter2d(double *x, int nx, int ny, int fsize_x, int fsize_y)
{
    /* quick return if kernel size is 1 */
    if(fsize_x <= 1 && fsize_y <= 1) return;

    /* make sure filter size is odd */
    if (!(fsize_x&1))
        fsize_x++;
    if (!(fsize_y&1))
        fsize_y++; 
    
    
    int i, j, k;

    int N = nx * ny;
    int total_fsize = fsize_x * fsize_y;

    int y_half = fsize_y / 2;
    int x_half = fsize_x / 2;
    
    
    /* select median function */ 
    /* need a good heuristic to select a median function  based on filter size */
    double (*med)(const double *x, const int N);
    if (total_fsize <= 9)
        med = med_small;
    else
    {
        //med = median_of_medians;
        med = median;
    }


    /* allocate temporary arrays */

    /* array to hold temp copy of x */
    double *y = malloc(sizeof(*y) * N);

    /* array of pointer to each row */
    double **rp = malloc(sizeof(**rp) * fsize_y);

    /* array to hold values for median */
    double *qq = malloc(sizeof(*qq) * total_fsize); 

    if (!y || !rp || !qq)
    {
        fprintf(stderr,"Error: median_filter2d: memory allocation failure\n");
        exit(1);
    }

    /* make a temp copy of x */
    memcpy(y, x, sizeof(*y) * N);

    /* loop through all rows */
    int nrows;
    int ncols;
    int rstart;
    int cstart;
    for (i=0;i<ny;i++)     
    {
        nrows = fsize_y;
        rstart = i-y_half;

        if (i<y_half)
        {
            nrows = (i+1) + y_half;
            rstart = 0;
        }

        if (i>=ny-y_half)
        {
            nrows = (ny - i) + y_half;
        }


        /* update row pointers */
        for (j=0;j<nrows;j++)
        {
            rp[j] = y + (j+rstart) * nx;
        }



        /* loop through each element in row */
        for (j=0;j<nx;j++)
        {
            ncols = fsize_x;
            cstart = j - x_half;
            if (j<x_half)
            {
                ncols = (j+1) + x_half;
                cstart = 0;
            }
            if (j>=nx-x_half)
            {
                ncols = (nx - j) + x_half;
            }

            

            for (k=0;k<nrows;k++)
            {
                /* copy data into temp array */
                memcpy(qq+k*ncols, rp[k]+cstart, sizeof(*qq) * ncols);
                
                /* increment row pointers */
                //++rp[k];
            }
            /* take median */

            
            x[j+ i*nx] = med(qq,nrows * ncols);
        }
    }




    /* free allocated memory */
    free(qq);
    free(rp);
    free(y);
}


void median_filter3d(double *x, int nx, int ny, int nz, int fsize_x, int fsize_y, int fsize_z)
{
    fprintf(stderr,"Warning: median_filter3d not yet implemented\n");
    fprintf(stderr,"Median filtering not performed\n");
    return;
}


//int main()
//{
//    //double x[] = {2,3,5,9,4,1,3,6,
//    //              4,2,5,6,1,8,1,0,
//    //              0,0,1,4,6,7,9,2,
//    //              9,8,1,3,4,6,0,8,
//    //              2,3,6,7,3,2,1,0 };
//    //int nx = 8;
//    //int ny = 5;
//    //int fsize = 5;
//
//    int fsize = 9;
//    int nx = 600;
//    int ny = 200;
//    int N = nx*ny;
//    double *x = malloc(sizeof(*x) * N);
//    int i;
//    for (i=0;i<N;i++)
//        x[i] = (double)rand() / (double)RAND_MAX;
//
//
//    
//
//    median_filter2d(x,nx,ny,fsize,fsize);
//
//    //int i,j;
//    //for (i=0;i<ny;i++)
//    //{
//    //    for (j=0;j<nx;j++)
//    //    {
//    //        printf("%.1f,",x[j+i*nx]);
//    //    }
//    //    printf("\n");
//    //}
//
//
//
//    return 0;
//}
//