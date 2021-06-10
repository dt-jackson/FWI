#include <stdlib.h>
#include <stdio.h>

#include <math.h>


void read_source(char *fname, int N, double *src)
{
    FILE *fid = fopen(fname,"r");
    if (!fid)
    {
        fprintf(stderr, "Error opening file %s\n",fname);
        exit(1);
    }

    char *line;
    size_t line_len = 0;;
    ssize_t char_read;

    int count = 0;
    double t;
    char *next;
    while ((char_read = getline(&line,&line_len,fid)) != -1)
    {
        t = strtod(line,&next);
        src[count] = strtod(next, NULL);
        count++; 
        line_len = 0;
    }

    free(line);
    fclose(fid);
}


void write_source(char *fname, double *src, int N, double t0, double dt)
{

    FILE *fid = fopen(fname,"w");
    if (!fid)
    {
        fprintf(stderr, "Error opening file %s\n",fname);
        exit(1);
    }

    int i;
    for (i=0;i<N;i++)
    {
        fprintf(fid, "%.16E    %.16E\n", t0+i*dt, src[i]);
    }

    fclose(fid);
}


void define_gaussian(double *g, int N, double t0, double dt, double tshift, double fmax)
{
    /* define a gaussian approximation to an impulse function */



    /* calculate nyquist frequency */
    double f_nyq = 1.0 / (2.0 * dt);

    /* determine fmax as fraction of nyquist */
    fmax = fmax / f_nyq;

    
    double sigma, sigma_f;

    /* limit for gaussian cutoff in freq domain */
    /* filter response at cutoff will be 1/c    */
   double c = 4.0;

    //sigma_f = fmax / sqrt(2.0 * log(c));
    sigma_f = fmax;

    sigma = 1.0 / sigma_f;

    sigma = 0.005;

    double t;
    double arg;
    double p1 = 1.0 / (sigma * sqrt(2.0 * M_PI));

    int i;
    for (i=0;i<N;i++)
    {
        t = (i*dt + t0);
        arg = (t - tshift) / sigma;

        g[i] = p1 * exp(- 0.5 * (arg * arg));
    }







}