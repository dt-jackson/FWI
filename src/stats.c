#include <stdlib.h>     /* malloc, qsort, exit      */
#include <stdio.h>      /* fprintf                  */
#include <string.h>     /* memcpy                   */
#include <math.h>       /* pow, sqrt, fabs          */

#include "stats.h"



int comp_double_ascending(const void *x1, const void *x2)
{
    /* for use with qsort           */
    /* puts smallest element first  */

    /* cast values as doubles */
    double a = *(double*)x1;
    double b = *(double*)x2;

    /* a is smaller than b, therefore a goes before b */
    if (a < b) return -1;

    /* a is larger than b, therefore a goes after b */
    if (a > b) return 1;

    /* otherwise, equal */
    return 0;
}

int comp_double_descending(const void *x1, const void *x2)
{
    /* for use with qsort           */
    /* puts largest element first   */

    /* cast values as double */
    double a = *(double*)x1;
    double b = *(double*)x2;
    
    /* a is larger than b, therefore a goes before b */
    if (a > b) return -1;


    /* a is smaller than b, therefore a goes after b */
    if (a < b) return 1;

    /* otherwise, equal */
    return 0;
}

int double_eq(const void * x1, const void *x2)
{
    /* returns 1 if two double are equal within STAT_TOL */
    /* returns 0 otherwise                               */

    /* cast values as double */
    double a = *(double*)x1;
    double b = *(double*)x2;

    /* check for equality */
    if (fabs(a - b) < STAT_TOL) return 1;

    /* return 0 if not equal */
    return 0;
}




double median(const double *x, const int N)
{
    /* calculate the median of the values in x */

    /* in order to avoid side effects, make a temp copy of input array */ 
    double *y = malloc(sizeof(*y) * N);
    if (!y)
    {
        fprintf(stderr, "Error: stats.c, median: memory allocation failure\n");
        exit(1);
    }
    memcpy(y, x, sizeof(*y) * N);
    

    double med;

    /* sort array */
    qsort(y, N, sizeof(*y), comp_double_ascending);

    /* determine if the number of elements is odd or even */
    int odd = N % 2; /* 0 if even, 1 if odd */

    int ind1, ind2;
    if (odd)
    {
        /* if N is odd, pick middle element */
        ind1 = 0.5 * (N-1);
        med = y[ind1];
    }
    else
    {
        /* if N is even, average two middle elements */
        ind1 = 0.5 * N;
        ind2 = ind1 -1;
        med = 0.5 * (y[ind1] + y[ind2]);
    }

    free(y);
    return med;
}


double mag_median(const double *x, const int N)
{
    /* calculates the median of the magnitudes of the values in x*/

    /* in order to avoid side effects, make a temp copy of input array */ 
    double *y = malloc(sizeof(*y) * N);
    if (!y)
    {
        fprintf(stderr, "Error: stats.c, median: memory allocation failure\n");
        exit(1);
    }

    int i;
    for (i=0;i<N;i++)
    {
        y[i] = fabs(x[i]);
    } 
    

    double med;

    /* sort array */
    qsort(y, N, sizeof(*y), comp_double_ascending);

    /* determine if the number of elements is odd or even */
    int odd = N % 2; /* 0 if even, 1 if odd */

    int ind1, ind2;
    if (odd)
    {
        /* if N is odd, pick middle element */
        ind1 = 0.5 * (N-1);
        med = y[ind1];
    }
    else
    {
        /* if N is even, average two middle elements */
        ind1 = 0.5 * N;
        ind2 = ind1 -1;
        med = 0.5 * (y[ind1] + y[ind2]);
    }

    free(y);
    return med;




}



double mean(const double *x, const int N)
{
    /* calculate the mean of the values in x */

    double m = 0.0;
    int i;

    /*sum elements */
    for (i=0;i<N;i++)
        m += x[i];
    
    /* divide by total number of elements and return */
    return m / N;
}





double percentile(const double *x, const int N, const int p)
{
    /* calculate the p-th percentile of the values in x */
    /* p should be between 0 and 100                    */

    /* in order to avoid side effects, make a temp copy of input array */ 
    double *y = malloc(sizeof(*y) * N);
    if (!y)
    {
        fprintf(stderr, "Error: stats.c, median: memory allocation failure\n");
        exit(1);
    }
    memcpy(y, x, sizeof(*y) * N);
    


    /* sort array */
    qsort(y, N, sizeof(*y), comp_double_ascending);

    /* get index using nearest rank method */
    int ind = (int)ceil(((double)p/100.0) * (double)N);

    /* subtract 1 because indexing starts at 0 */
    double m = y[ind-1];

    free(y);
    return m;

}





void insertion_sort(double *x, int N)
{
    /* use this when N is small */

    int  i, j;
    double t;
    for (i=1;i<N;i++)
    {
        t = x[i];
        for (j=i; j>0 && x[j-1] > t;j--)
            x[j] = x[j-1];
        x[j] = t;
    }
}


void isort_alt(double *x, int N)
{
    /* alternate insertion sort */
    int i, j;
    double t;
    for (i=1; i<N;i++)
    {
        for (j=i; j<0 && x[j-1] > x[j]; j--)
        {
            t = x[j];
            x[j] = x[j-1];
            x[j-1] = t;
        }
    }
    
}


double median_of_medians_n(const double *x, const int N, const int j)
{
    /* not quite working properly */
    fprintf(stderr,"median_of_medians still has a bug\n");
    fprintf(stderr,"Use median or med_small instead\n");
    return median(x,N);

    /* for large N, this is much more efficient */


    double *y = malloc(sizeof(*y) * N);
    memcpy(y, x, sizeof(*y) * N);

    int sublist_size = 5;
    int med_index = 2;
    int no_sublists = (N / sublist_size) + 1;
    int last_list_sz = sublist_size - ((sublist_size * no_sublists) - N);
    int odd = (last_list_sz & 1) - 1; /* 0 if odd, -1 if even */
    int last_list_med_index = last_list_sz/2;


    double *med = malloc(sizeof(*med) * no_sublists);
    double pivot;

    int i;
    double *list = y;
    for (i=0;i<no_sublists-1;i++)
    {
        /* pick medians of sublists */
        insertion_sort(list, sublist_size);
        med[i] = list[med_index];
        list += sublist_size;
    }
    /* do last list */
    if (last_list_sz)
    {
        insertion_sort(list, last_list_sz);
        med[i] = list[last_list_med_index];
        /* an extra add and mutliply here if last_list_sz is odd, but  */
        /* we avoid branching - need to see which is faster            */
        med[i] += list[last_list_med_index + odd];
        med[i] *= 0.5;
    }

    if (no_sublists < sublist_size)
    {
        insertion_sort(med, no_sublists);
        pivot = med[no_sublists/2];
    }
    else 
    {
        pivot = median_of_medians_n(med, no_sublists, no_sublists/2);
    }

    /* partitionling step */
    double *low = malloc(sizeof(*low) * N);     /* values in y less than pivot*/
    double *high = malloc(sizeof(*high) * N);   /* values in y greater than pivot */
    int len_low = 0;
    int len_hi = 0;
    for (i=0;i<N;i++)
    {
        if (y[i] < pivot)
            low[len_low++] = y[i];
        else if (y[i] > pivot)
            high[len_hi++] = y[i];
    }


    int k = len_low; /* length of low */
    if (j < k)
    {
        return median_of_medians_n(low, len_low, j);
    } 
    else if (j > k)
    {
        return median_of_medians_n(high, len_hi, j-k-1);
    }
    else 
    {
        return pivot;
    }

    free(low);
    free(high);
    free(med);
    free(y);
    return pivot;
}

double median_of_medians(const double *x, const int N)
{
    return median_of_medians_n(x,N, N/2);
}



double med_small(const double *x, const int N)
{
    /* for when N<=~10  */
    double y[25];
    memcpy(y, x, sizeof(*x) * N);
    insertion_sort(y,N);
    int odd = N & 1;
    if (odd)
        return y[N/2];
    else 
        return 0.5 * (y[N/2] + y[N/2 -1]);
}


