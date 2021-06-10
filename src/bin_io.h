#ifndef BIN_IO_H
#define BIN_IO_H

#include <stdio.h>


FILE * open_bin_read(char *fname);
FILE * open_bin_write(char *fname);

int get_num_records(FILE *binfile);

int read_bin_float(FILE *binfile,float *values);
void write_bin_float(int N,float *values,FILE *binfile);

int read_NSPEC_IBOOL(FILE *binfile,int *ibool);




#endif