
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <color.h>

void v_compute_min_max(unsigned char* v_in, long size, long offset, unsigned char* v_min, unsigned char* v_max)
{
	int i;
	unsigned char v, vmin, vmax;

	vmin = vmax = v_in[offset];

	for (i = offset+4 ; i < size*4; i+=4)
	{
		v = v_in[i];
		if (v < vmin) {
			vmin = v;
		}
		else if (v > vmax) {
			vmax = v;
		}
	}
	// Format compact : R-G-B-A-R-G-B-A..

	/*for (i = offset + 4; i < size * 4; i += 4)
	{
		v = v_in[i*4+offset];
		if (v < vmin) {
			vmin = v;
		}
		else if (v > vmax) {
			vmax = v;
		}
	}*/

	*v_min = vmin;
	*v_max = vmax;
}

typedef struct _AREA_216_ 
{
	unsigned char bary_r, bary_g, bary_b;
	long nb_points;
	long sum_r, sum_g, sum_b;
} AREA_216;



static int v_to_no(double v, double v_min, double v_max, int nb)
{
	int no;
	if (v_max == v_min) return 0;
	no = (int)floor((v-v_min)/((v_max - v_min) / (double)nb));
	if (no < 0) no = 0; else if (no >= nb) no = nb-1;
	return no;
}



void img_reduce_colors_to_216(unsigned char * img_color_in,
															long width, long height,
															unsigned char * img_color_out)
{
	// Max et min RGB
	unsigned char min_rgb[3], max_rgb[3], * img_no = NULL;
	
	int k, size = width * height, i, no_rgb[3], no;
	AREA_216* area216 = NULL;
	
	
	/*
	unsigned char* ptr1 = NULL, * ptr2 = NULL, * ptr3 = NULL;
	ptr1 = &(min_rgb[2]);
	ptr2 = min_rgb + 2;
	ptr3 = &min_rgb + 2;
	// ptr1 et ptr2 sont equivalentes
	*/
	
	for(k=0; k<3; k++)
		v_compute_min_max(img_color_in, size, k, &(min_rgb[k]), &(max_rgb[k]));
	area216 = (AREA_216*)calloc(216, sizeof(AREA_216));

	img_no = (unsigned char *)calloc(size, sizeof(unsigned char));
	for (i = 0; i < size; i++)
	{
		for (k = 0; k < 3 ; k++)
		{
			no_rgb[k] = v_to_no((double)img_color_in[i*4+k], (double)min_rgb[k], (double)max_rgb[k], 6);
		}
		no = no_rgb[0] + 6 * no_rgb[1] + 36 * no_rgb[2];
		img_no[i] = (unsigned char)no;
		area216[no].nb_points++;
		area216[no].sum_r += (long)img_color_in[i * 4];
		area216[no].sum_g += (long)img_color_in[i * 4 + 1];
		area216[no].sum_b += (long)img_color_in[i * 4 + 2];
	}

	for (i = 0; i < 216; i++)
	{
		if (area216[i].nb_points)
		{
			area216[i].bary_r = (unsigned char)floor((double)area216[i].sum_r / (double)area216[i].nb_points + 0.5);
			area216[i].bary_g = (unsigned char)floor((double)area216[i].sum_g / (double)area216[i].nb_points + 0.5);
			area216[i].bary_b = (unsigned char)floor((double)area216[i].sum_b / (double)area216[i].nb_points + 0.5);
		}
	}

	for (i = 0; i < size; i++)
	{
		no = img_no[i];
		img_color_out[i * 4] = area216[no].bary_r;
		img_color_out[i * 4 + 1] = area216[no].bary_g;
		img_color_out[i * 4 + 2] = area216[no].bary_b;
	}



	free(area216);
}




