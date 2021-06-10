#include <stdlib.h>
#include <stdio.h>

#include "read_parameters.h"

#include "seis_data.h"


void read_src_loc(parameters *p, seis_data *d)
{
    int num_src = d->num_src;
    
    /* open src loc file */ 
    FILE *src_file = fopen(p->src_loc_file, "r");
    if (!src_file)
    {
        fprintf(stderr, "Erroring opening file %s\n",p->src_loc_file);
        exit(2);
    }

    char *line;
    size_t line_len = 0;
    int k = 0;

    while (getline(&line,&line_len,src_file) != -1)
    {
        remove_comments(line);
        remove_white_space(line);
        if (line[0] != '\0')
        {
            d->src_loc[k] = atof(line);
            k++;
        }
    }

    printf("%d source locations read\n",k);
    if (k != num_src)
    {
        fprintf(stderr, "Source locations not read properly\n");
        exit(2);
    }

    free(line);
    fclose(src_file);    
}


void read_rec_loc(parameters *p, seis_data *d)
{
    int num_rec = d->num_rec;
    
    /* open src loc file */ 
    FILE *rec_file = fopen(p->rec_loc_file, "r");
    if (!rec_file)
    {
        fprintf(stderr, "Erroring opening file %s\n",p->rec_loc_file);
        exit(2);
    }

    char *line;
    size_t line_len = 0;
    int k = 0;

    while (getline(&line,&line_len,rec_file) != -1)
    {
        remove_comments(line);
        remove_white_space(line);
        if (line[0] != '\0')
        {
            d->rec_loc[k] = atof(line);
            k++;
        }
    }

    printf("%d receiver locations read\n",k);
    if (k != num_rec)
    {
        fprintf(stderr, "Receiver locations not read properly\n");
        exit(2);
    }
    free(line);
    fclose(rec_file);    
}