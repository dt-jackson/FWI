#ifndef READ_SU_H
#define READ_SU_H

#include <stdio.h>

#include "seis_data.h"


void read_su_data_only(char *file_name, seis_data *d, int src_no);

/* function read_su is not quite complete yet - doesn't fully read headers */
void read_su(char *file_name,seis_data *d, int src_no, int *nstep);



void get_num_rec_su(char *file_name, int* num_samp, int *num_rec);



int read_1byte_int(FILE *fn);
int read_2byte_int(FILE *fn);
int read_4byte_int(FILE *fn);
double read_fp_single(FILE *fn);
double read_fp_double(FILE *fn);

#endif