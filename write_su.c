#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "bin_io.h"
#include "seis_data.h"
#include "write_su.h"

void write_su_data_only(char *file_name, seis_data *d, int src_no)
{
    long i, j;
    int header_sz = 240; /* size of header in bytes     */
    int float_sz = 4;    /* size of floating point data */

    long ns = d->num_samples; /* number of samples */

    /* open file for writing */
    FILE *su_file = open_bin_write(file_name);

    /* copy data into temp single precision array */
    long tot_samp = d->num_samples * d->num_rec;
    float *temp = malloc(sizeof(*temp) * tot_samp);
    if (!temp)
    {
        fprintf(stderr,"Error: write_su_data_only: memory allocation failure\n");
        exit(1);
    }
    
    for (i=0;i<tot_samp;i++)
    {
        if (src_no < 0)
            temp[i] = 0.0;
        else
            temp[i] = d->traces[i + src_no*tot_samp];
    }

    /* write data to file, skipping over headers */
    fseek(su_file,0,SEEK_SET);
    for (i=0;i<d->num_rec;i++)
    {
        fseek(su_file,header_sz,SEEK_CUR);
        fwrite(temp + i*ns,float_sz,ns,su_file);
    }

    /* free local array */
    free(temp);

    /* close file for writing */
    fclose(su_file);
}

void write_su(char *file_name, seis_data *d, int src_no)
{
    /* this is not entirely correct, but its good enough 
        for specfem2d */


    long i, j;
    int header_sz = 240; /* size of header in bytes     */
    int float_sz = 4;    /* size of floating point data */

    long ns = d->num_samples; /* number of samples */

    /* open file for writing */
    FILE *su_file = fopen(file_name,"wb");

    /* copy data into temp single precision array */
    long tot_samp = d->num_samples * d->num_rec;
    float *temp = malloc(sizeof(*temp) * tot_samp);
    
    for (i=0;i<tot_samp;i++)
    {
        if (src_no < 0)
            temp[i] = 0.0;
        else
            temp[i] = d->traces[i + src_no*tot_samp];
    }








    /* make sure we're at the beginning of the file */
    fseek(su_file,0,SEEK_SET);
    for (i=0;i<d->num_rec;i++)
    {
        /* write header - 240 bytes*/
        int32_t tracl = (int32_t)(i+1);
        int32_t offset = (int32_t)(d->rec_loc[i] - d->src_loc[src_no]);
        int32_t sx = (int32_t)d->src_loc[src_no];
        int32_t gx = (int32_t)d->rec_loc[i];
        int16_t ns_ = (int16_t)ns;
        int16_t dt_ = (int16_t)(d->dt * 1e6);


        /* write tracl at start */
        fwrite(&tracl,4,1,su_file);

        /* skip 32 bytes to get to offset */
        fseek(su_file,32,SEEK_CUR);
        fwrite(&offset,4,1,su_file);

        /* skip 32 bytes to get to sx */
        fseek(su_file,32,SEEK_CUR);
        fwrite(&sx,4,1,su_file);

        /* skip 4 bytes to get to gx */
        fseek(su_file,4,SEEK_CUR);
        fwrite(&gx,4,1,su_file);

        /* skip 30 bytes to get to ns */
        fseek(su_file,30,SEEK_CUR);
        fwrite(&ns_,2,1,su_file);

        /* skip 0 bytes to get to dt */
        fwrite(&dt_,2,1,su_file);

        /* skip 122 bytes to get to end of header */
        fseek(su_file,122,SEEK_CUR);



        /* write data */
        fwrite(temp + i*ns,float_sz,ns,su_file);
    }

    /* free local array */
    free(temp);

    /* close file for writing */
    fclose(su_file);



}