#include <stdlib.h>
#include <stdio.h>


int read_bin_float(FILE *binfile,float *values)
{
    /* move to the beginning of the file */
    fseek(binfile,0,SEEK_SET);

    /* Read opening tag - 4-byte integer */
    int otag; 
    fread((int*)&otag,4,1,binfile);
    
    /* caluclate numebr of values. Each value is a 4-byte float */
    int N = otag / 4;
    

    /* Read records - N 4-byte floating point values */
    fread((float*)values,4,N,binfile);

    /* Read closing tag - 4-byte integer */
    int ctag;
    fread((int*)&ctag,4,1,binfile);
    
    /* The closing tag should be the same as the opening tag */
    if (ctag != otag)
    {
        /* TODO: add error message here and exit*/
    }

    /* return number of values successfully read */
    return N;
}

void write_bin_float(int N,float *values,FILE *binfile)
{
    /* calculate number of bytes. each value is a 4 byte float */
    int num_bytes = N * 4;

    /* write opening tag - 4-byte integer */
    fwrite(&num_bytes,4,1,binfile);

    /* write data - N 4-byte floats */
    fwrite(values,4,N,binfile);
    
    /* write closing tag - 4-byte integer */
    fwrite(&num_bytes,4,1,binfile);
}

/*****************************************************************************/

int read_NSPEC_IBOOL(FILE *binfile,int *ibool) 
{
    /* Read NSPEC otag */
    int otag;
    fread(&otag,4,1,binfile);

    /* Read NSPEC */
    int NSPEC;
    fread(&NSPEC,4,1,binfile);

    /* Read NSPEC ctag */
    int ctag;
    fread(&ctag,4,1,binfile);

    /* verify that otag and ctag are the same */
    if (otag != ctag) 
    {
        /* TODO: add error message here and exit*/
    }

    /* Read ibool otag */
    fread(&otag,4,1,binfile);

    /* caluclate numebr of values. Each value is a 4-byte int */
    int N = otag / 4;

    /* Read ibool */
    fread(ibool,4,N,binfile);

    /* read ibool ctag */
    fread(&ctag,4,1,binfile);


    /* verify that otag and ctag are the same */
    if (otag != ctag) 
    {
        /* TODO: add error message here and exit*/
    }

    return NSPEC;
}


/*****************************************************************************/
/* Functions for getting the number of records from a bin file               */

int get_num_records(FILE *binfile) 
{
    /* This function gets the number of records from a single file */

    /* move to the beginnning of the file */
    fseek(binfile,0,SEEK_SET);

    /* read number of bytes in records - 4 byte integer */
    int N;
    fread((int*)&N,4,1,binfile); 

    /* divide by 4 because each value is a 4-byte float */
    return N/4;
}






/*****************************************************************************/
/* wrapper functions for opening bin files. This checks to make sure they    */
/* were properly opened                                                      */
FILE * open_bin_read(char *fname) 
{
    FILE *binfile = fopen(fname,"rb");
    if (binfile == NULL)
    {
        fprintf(stderr,"Error opening file %s\n",fname);
        exit(1);
    }
    return binfile;
}

FILE * open_bin_write(char *fname) 
{
    FILE *binfile = fopen(fname,"wb");
    if (binfile == NULL)
    {
        fprintf(stderr,"Error opening file %s\n",fname);
        exit(1);
    }
    return binfile;
}