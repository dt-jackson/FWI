#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "bin_io.h"
#include "seis_data.h"
#include "read_su.h"

/*****************************************************************************/

void read_su_data_only(char *file_name, seis_data *d, int src_no)
{

    /* header size in bytes */
	int head_sz = 240; 
    /* bytes per floating point sample */
    int float_sz = 4;

    long ns = d->num_samples; 
    long num_rec = d->num_rec;
    long ns_tot = num_rec * ns;

    /* open file */
    FILE *su_file = open_bin_read(file_name);

    /* make a temp single precision array to read the data in to */
    float *temp = malloc(sizeof(*temp)*ns_tot);
    if (!temp)
    {
        fprintf(stderr,"Memory allocation failure in read_su_data_only\n");
        exit(1);
    }

    /* start at the beginning of the file and read the data, */
    /* skipping over the headers                             */
    int i;
    fseek(su_file,0,SEEK_SET);
    for (i=0; i<num_rec; i++)
    {
        fseek(su_file,head_sz,SEEK_CUR);
        fread(temp+i*ns,4,ns,su_file);
    }

    /* copy data to double array in seis_data struct */
    for (i=0;i<ns_tot;i++)
    {
        d->traces[i+src_no*ns_tot] = temp[i];
    }
    
    /* free local array */
    free(temp);

    /* close file */
    fclose(su_file);
}

/*****************************************************************************/

/* this function is not quite complete yet - doesn't fully read headers */
void read_su(char *file_name,seis_data *d, int src_no, int *nstep)
{
    /* if the number of samples should be read from the file, set     */
    /* nstep to NULL                                                  */
    /* The number of samples is recorded in the file as a 2 byte int, */
    /* so it maxes out at 65535. However, the file can hold more than */
    /* that. In these cases, the number of samples should be passes   */
    /* in via nstep.                                                  */ 


    /* header size in bytes */
	int head_sz = 240; 

    /* number of bytes to skip to get number of samples */
	int ns_loc = 114;

    /* bytes per floating point sample */
    int float_sz = 4;

    /* open file */
    FILE *su_file = fopen(file_name,"rb");

    /* get file length */
    fseek(su_file,0,SEEK_SET);
    long beg = ftell(su_file);
    fseek(su_file,0,SEEK_END);
    long end = ftell(su_file);
    long flength  = end - beg;



    /* read number of samples */
    fseek(su_file,ns_loc,SEEK_SET);
    int ns = read_2byte_int(su_file);
    if (nstep != NULL) ns = *nstep;

    /* read the time step - comes right after num samples */
    double dt = (double)read_2byte_int(su_file); /* this is in micro sec */
    dt *= (1.0e-6); /* convert from usec to sec */

    /* calculate number of traces */
    int num_rec = flength / (head_sz + float_sz * ns);

    /* make a temp single precision array to read the data in to */
    float *temp = malloc(sizeof(*temp)*ns*num_rec);
    if (!temp)
    {
        fprintf(stderr,"Memory allocation failure int read_su\n");
        exit(1);
    }


    /* go back to the beginning and read the data, skipping over the headers */
    int i;
    fseek(su_file,0,SEEK_SET);
    for (i=0;i<num_rec;i++)
    {
        fseek(su_file,head_sz,SEEK_CUR);
        fread(temp+i*ns,4,ns,su_file);
    }

    /* copy data to double array in seis_data struct */
    int ns_tot = num_rec * ns;
    for (i=0;i<ns_tot;i++)
    {
        d->traces[i+src_no*ns_tot] = temp[i];
    }
    

    free(temp);
    /* close file */
    fclose(su_file);
}

/*****************************************************************************/

void get_num_rec_su(char *file_name, int* num_samp, int *num_rec)
{
    int header_len = 240;   /* length of each header in bytes   */
    int ns_offset = 114;    /* offset to number of samples      */
    int float_sz = 4;       /* bytes per floating point sample  */

    /* open file */
    FILE *su_file = fopen(file_name,"rb");

    /* get number of samples */
    fseek(su_file,ns_offset,SEEK_SET);
    *num_samp = read_2byte_int(su_file);
    
    /* get file length */
    fseek(su_file,0,SEEK_SET);
    long beg = ftell(su_file);
    fseek(su_file,0,SEEK_END);
    long end = ftell(su_file);
    long flen = end - beg;

    *num_rec = flen / (header_len + float_sz * (*num_samp));

    /* close file */
    fclose(su_file);
    
    return;
}




/*****************************************************************************/

double read_fp_single(FILE *fn)
{
    float value;
    fread((float*)&value,4,1,fn);
    return (double)value;
}

double read_fp_double(FILE *fn)
{
    double value;
    fread((double*)&value,8,1,fn);
    return value;
}

int read_1byte_int(FILE *fn)
{
    int8_t value;
    fread((int8_t*)&value,1,1,fn);
    return value;
}

int read_2byte_int(FILE *fn)
{
    int16_t value;
    fread((int16_t*)&value,2,1,fn); 
    return (int)value;
}

int read_4byte_int(FILE *fn)
{
    int32_t value;
    fread((int32_t*)&value,4,1,fn); 
    return (int)value;
}