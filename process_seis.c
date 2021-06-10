#include <stdlib.h>
#include <stdio.h>
#include <math.h>



#include "seis_data.h"


void sr_dist_mute(seis_data *s, double sr_min, double sr_max)
{
    /* mute traces closer than sr_min or farther than sr_max */

    double *trace_ptr;
    double sr_dist;
    int num_src = s->num_src;
    int num_rec = s->num_rec;
    int nstep = s->num_samples;
    int rec_step = num_rec * nstep;
    int i, j, k;

    /* loop through all sources */
    for (i=0; i<num_src; i++)
    {
        /* loop through all receivers */
        for (j=0; j<num_rec; j++)
        {
            /* pointer to start of this  trace */
            trace_ptr = s->traces + i*rec_step + j*nstep;

            /* calculate source-receiver distance */
            sr_dist = fabs(s->src_loc[i] - s->rec_loc[j]);

            /* if sr_dist is outside of bounds, set trace to zero */
            if (sr_dist < sr_min || sr_dist > sr_max)
            {
                for (k=0;k<nstep;k++)
                {
                    trace_ptr[k] = 0.0;
                }
            }
        }

    }
}