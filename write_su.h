#ifndef WRITE_SU_H
#define WRITE_SU_H

#include "seis_data.h"

void write_su_data_only(char *file_name, seis_data *d, int src_no);
void write_su(char *file_name, seis_data *d, int src_no);

#endif