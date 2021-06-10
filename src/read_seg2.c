#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>

#include "seis_data.h"
#include "read_seg2.h"


void read_seg2(const char *file_name, seis_data *dd)
{

    /* open seg2 file */
    FILE *seg2 = fopen(file_name,"rb");
    if (!seg2)
    {
        fprintf(stderr,"Error opening file %s\n",file_name);
        exit(-1);
    }


	/* read file descriptor block - 2 bytes */
	int16_t fdb;
	fread(&fdb,2,1,seg2);


	/* read revision number - 2 bytes */
	int16_t rev_no;
	fread(&rev_no, 2,1,seg2);

	/* size of trace pointer sub-block - 2 bytes */
	int16_t tpsb;
	fread(&tpsb, 2,1,seg2);
    int M = tpsb;


	/* number of traces - 2 bytes */
	int16_t num_rec;
	fread(&num_rec, 2,1,seg2);
    int N = num_rec;

	/* size of string terminator - 1 byte */
	char st_sz;
	fread((char*)&st_sz, 1,1,seg2);

	/* first string terminator char - 1 byte */
	char st1;
	fread((char*)&st1, 1,1,seg2);

    /* second strings terminator - 1 byte */
	char st2;
	fread((char*)&st2, 1,1,seg2);

	/* size of line terminator - 1 byte */
	char lt_sz;
	fread((char*)&lt_sz, 1,1,seg2);

	/* first line terminator character - 1 byte */
	char lt1;
	fread((char*)&lt1, 1,1,seg2);

	/* second line terminator character -             */
    /* empty if size of line terminator is 1 - 1 byte */
	char lt2;
	fread((char*)&lt2, 1,1,seg2);

	/* 18 reserved bytes */
    fseek(seg2,18,SEEK_CUR);

	//trace pointers 
    int32_t *trace_pointers = malloc(sizeof(*trace_pointers) * N);
	if (!trace_pointers)
	{
		fprintf(stderr,"Error: read_seg2: memory allocation failure\n");
		exit(1);
	}
    fread((int32_t*)trace_pointers,4,N,seg2);

	/* bytes left in trace pointer sub block */
	int bl = M - 4 * N;
    fseek(seg2,bl,SEEK_CUR);

    /* read file strings */
	seg2_header tmp_h;
	int s_offset = read_strings(seg2,&tmp_h);


	/* skips */
	int nskips = trace_pointers[0] - (32 + M + s_offset + 2);
	fseek(seg2, nskips, SEEK_CUR);


	/* read traces */
	int ii;
	void *tmp_data;
	for (ii=0;ii<N;ii++)
	{
		/* read trace ID - 2 bytes*/
		int16_t traceID;
		fread(&traceID, 2, 1, seg2);

		/* read block size - 2 bytes */
		int16_t blockSize;
		fread(&blockSize, 2, 1, seg2);

		/* read size of data block - 4 bytes */
		int32_t dataSize;
		fread(&dataSize, 4, 1, seg2);

		//read number of samples - 4 bytes
		int32_t num_samp;
		fread(&num_samp, 4, 1, seg2);

		//read data format code - 1 byte
		int8_t data_format_code;
		fread(&data_format_code, 1, 1, seg2);

		//19 bytes reserved
		fseek(seg2, 19, SEEK_CUR);

		/* for first trace, allocate space for data */
		if (ii == 0)
		{
			/* allocate temporary space */
			if (data_format_code == 4)
			{
				/* 4 byte floating point */
				tmp_data = malloc(sizeof(float) * num_samp);
			}
			else if (data_format_code == 5)
			{
				/* 8 byte floating point */
				tmp_data = malloc(sizeof(double) * num_samp);
			}
			else 
			{
				/* unknown format code */
				fprintf(stderr, "Error: read_seg2: Unknow data format code (%d) for file %s\n",(int)data_format_code, file_name);
				exit(1);
			}


			int num_src = 1;
			allocate_seis_data(dd, num_src, num_rec, num_samp);
		}

		/* read strings */
		seg2_header th;
		int tlen = read_strings(seg2, &th);

		/* copy data over */
		dd->rec_loc[ii] = th.rec_loc;
		if (ii == 0)
		{
			dd->dt = th.dt;
			dd->delay = th.delay;
			dd->src_loc[ii] = th.src_loc;
		}

		/* skips */
		nskips = blockSize - (32 + tlen + 2);
		fseek(seg2, nskips, SEEK_CUR);

		/* read data into temp array, apply descaling factor, and copy to output struct*/
		int j;
		if (data_format_code == 4)
		{
			/* single precision */
			fread(tmp_data, sizeof(float), num_samp, seg2);

			for (j=0;j<num_samp;j++)
			{
				/* multiply by 0.001 to convert from mV to Volts */
				dd->traces[ii*num_samp + j] = 0.001 * th.descaling_factor * ((float*)tmp_data)[j];
			}

		}
		else if (data_format_code == 5)
		{
			/* double precision */
			fread(tmp_data, sizeof(double), num_samp, seg2);
			
			for (j=0;j<num_samp;j++)
			{
				/* multiply by 0.001 to convert from mV to Volts */
				dd->traces[ii*num_samp + j] = 0.001 * th.descaling_factor * ((double*)tmp_data)[j];
			}
		}
		else
		{
			/* unknown format code */
			fprintf(stderr, "Error: read_seg2: Unknow data format code (%d) for file %s\n",(int)data_format_code, file_name);
			exit(1);
		}
	}


	
	free(tmp_data);
	free(trace_pointers);

    /* close seg2 file */
    fclose(seg2);
}






#define MAX_STRING 4096
int read_strings(FILE *seg2, seg2_header *h)
{
	/* read the strings and match a limited number of keys */

	/* all these keys are doubles */
	int no_keys = 5;
	const char *keys[no_keys];
	keys[0] = "DELAY";
	keys[1] = "DESCALING_FACTOR";
	keys[2] = "RECEIVER_LOCATION";
	keys[3] = "SOURCE_LOCATION";
	keys[4] = "SAMPLE_INTERVAL";

	double values[no_keys];

	int i, j;
	for (i=0;i<no_keys;i++)
		values[i] = 0.0;

	/* read offset - 2 bytes */
	int16_t s_offset;
	fread(&s_offset,2,1,seg2);
	int total_offset = (int)s_offset;

	int count = 0;
	while (s_offset != 0)
	{
		/* read string */
		char file_s[MAX_STRING];
		fread(file_s, 1, s_offset - 2, seg2);

		/* make sure string is upper case */
		i = 0;
		while (file_s[i])
		{
			file_s[i] = toupper(file_s[i]);
			i++;
		}

		/* split string at space - only do this once*/
		char *kk, *vv;
		kk = strtok(file_s,	" \t");
		vv = strtok(NULL,   " \t");

		/* compare string against keys */
		for (i=0;i<no_keys;i++)
		{
			if((strcmp(kk,keys[i]) == 0) && vv)
			{
				/* key match and value is not NULL */
				/* convert to double */
				values[i] = atof(vv);
			}
		}

		/* read the next offset */
		fread(&s_offset,2,1,seg2);
		total_offset += (int)s_offset;

	}

	/* copy data to outputs */
	h->delay 			= values[0];
	h->descaling_factor = values[1];
	h->rec_loc 			= values[2];
	h->src_loc 			= values[3];
	h->dt		 		= values[4];



    return total_offset;
}



//int main(int argc, char **argv)
//{
//	seis_data dd;
//
//	/* read data */
//	read_seg2(argv[1], &dd);
//
//	/* print header */
//	printf("dt: %f\n", dd.dt);
//	printf("num_samp: %d\n", dd.num_samples);
//	printf("num_src: %d\n", dd.num_src);
//	printf("num_rec: %d\n", dd.num_rec);
//	int i, j;
//	for (i=0;i<dd.num_src;i++)
//		printf("src_loc %d: %f\n",i+1, dd.src_loc[i]);	
//	printf("\n");
//	for (i=0;i<dd.num_rec;i++)
//		printf("rec_loc %d: %f\n",i+1, dd.rec_loc[i]);	
//	printf("\n");
//
//	/* free data */
//	free_seis_data(&dd);
//	return 0;
//}