#ifndef STATS_HEADER
#define STATS_HEADER

/* tolerance for double equality comparisons */
#define STAT_TOL 1e-6

/* comparisons of two doubles - use with qsort and bsearch */
int comp_double_ascending(const void *x1, const void *x2);
int comp_double_descending(const void *x1, const void *x2);

/* basic univariate stat function */
double median(const double *x, const int N);
double mag_median(const double *x, const int N);
double mean(const double *x, const int N);
double percentile(const double *x, const int N, const int p);


void insertion_sort(double *x, int N);

/* for large N, this will be faster (although it's not optimized yet) */
/* median_of_medians_n can return any percentile.                     */
/* this doesn't quite  work yet */
double median_of_medians_n(const double *x, const int N, const int j);

/* special case of median_of_medians_n where j = N/2                */
double median_of_medians(const double *x, const int N);

/* this may be faster when N is small (less than 10 -15) and odd */
double med_small(const double *x, const int N);

#endif