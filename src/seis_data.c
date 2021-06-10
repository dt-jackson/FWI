#include <stdlib.h>
#include <stdio.h>


#define NDEBUG
#ifdef DEBUG
#undef NDEBUG
#endif
#include <assert.h>


#include "seis_data.h"



void allocate_seis_data(seis_data *t, int num_src, int num_rec, int num_samp) 
{
    assert(t != NULL);

    t->num_src = num_src;
    t->num_samples = num_samp;
    t->num_rec = num_rec;

    t->src_loc = malloc(sizeof(*t->src_loc)*num_src);
    t->rec_loc = malloc(sizeof(*t->rec_loc)*num_rec);
    t->traces  = malloc(sizeof(*t->traces)*num_src*num_rec*num_samp);




    if (!t->src_loc || !t->rec_loc || !t->traces)
    {
        fprintf(stderr,"Memory allocation failure in allocate_seis_data\n");
        exit(1);
    }



}

void free_seis_data(seis_data *t)
{
    if (!t)
        return;

    if (t->src_loc)
    {
        free(t->src_loc);
        t->src_loc = NULL;
    }

    if (t->rec_loc) 
    {
        free(t->rec_loc);
        t->rec_loc = NULL;
    }

    if (t->traces)
    {
        free(t->traces);
        t->traces = NULL;
    }
}


void copy_seis_data(seis_data *dest, seis_data *src)
{
    int num_src = src->num_src;
    int num_rec = src->num_rec;
    int nstep = src->num_samples;
    int total_step = num_src * num_rec * nstep;
    free_seis_data(dest);
    allocate_seis_data(dest,num_src,num_rec,nstep);

    int i;
    
    /* copy source locations */
    for (i=0;i<num_src;i++)
        dest->src_loc[i] = src->src_loc[i];

    /* copy receiver locations */
    for (i=0;i<num_rec;i++)
        dest->rec_loc[i] = dest->rec_loc[i];

    /* copy data */
    for (i=0;i<total_step;i++)
        dest->traces[i] = src->traces[i]; 


}


void print_seis(FILE *ofile, seis_data *d)
{
    /* print information */
    fprintf(ofile,"Number of sorces: %d\n",d->num_src);
    fprintf(ofile,"Number of receivers: %d\n",d->num_rec);
    fprintf(ofile,"Number of steps: %d\n",d->num_samples);
    fprintf(ofile,"Time step: %f\n",d->dt);

    /* print data */
    int num_src = d->num_src;
    int num_rec = d->num_rec;
    int nstep = d->num_samples;
    int i, j, k;
    for (i=0;i<num_src;i++)
    {
        fprintf(ofile,"Source Number %d\n",i+1);
        for (j=0;j<nstep;j++)
        {
            fprintf(ofile,"%d,",j);
            for (k=0;k<num_rec;k++)
            {
                fprintf(ofile,"%E,",d->traces[i*nstep*num_rec + k*nstep + j]);
            }
            fprintf(ofile,"\n");
        }
    }



}