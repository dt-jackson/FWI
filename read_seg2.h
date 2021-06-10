#ifndef READ_SEG2_H
#define READ_SEG2_H

#include <stdio.h>
#include "seis_data.h"

typedef struct seg2_header
{
    double delay;
    double descaling_factor;
    double src_loc;
    double rec_loc;
    double dt;

} seg2_header;


void read_seg2(const char *file_name, seis_data *dd);
int read_strings(FILE *seg2, seg2_header *h);

#endif