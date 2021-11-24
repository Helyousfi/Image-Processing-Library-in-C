

#include <rotate.h>


typedef struct _RGBA_
{
	unsigned char r;
	unsigned char g;
	unsigned char b;
	unsigned char a;
} RGBA;

void img_rotate(void * img_in, long width, long height,
								long type_data,
								void * img_out)
{
	/*
	int x,y;
	double* img_64f_in = (double*)img_in, * img_64f_out = (double*)img_out;
	RGBA* img_rgba_in = (RGBA*)img_in;
	RGBA* img_rgba_out = (RGBA*)img_out;
	if (type_data == 0)
	{
		for (y = 0; y < height; y++)
		{
			for (x = 0; x < width; x++)
			{
				//((double *)img_out)[(height - 1) - y + x * height] = ((double*)img_in)[x + y * width];
				img_64f_out[x * height + height - 1 - y] = img_64f_in[x + y * width];
				
			}
		}
	}

	else if (type_data == 1)
	{
		for (y = 0; y < height; y++)
		{
			for (x = 0; x < width; x++)
			{
				//((double *)img_out)[(height - 1) - y + x * height] = ((double*)img_in)[x + y * width];
				img_rgba_out[x * height + height - 1 - y] = img_rgba_in[x + y * width];
			}
	}

	*/
	
	int i, j, k, t;
	unsigned char* img_uchar_in = (unsigned char*)img_in;
	unsigned char* img_uchar_out = (unsigned char*)img_out;;
	if (type_data == 0) t = 8; else if (type_data == 1) t = 4;

	for (j = 0; j < height; j++)
	{
		for (i = 0; i < width; i++)
		{
			for (int k = 0; k < t; k++)
			{
				img_uchar_out[(i * height + height - 1 - j) * t + k] = img_uchar_in[(j * width + i) * t + k]; //(j * width + i) * t + k = j*width*t + i*t + k
			}
		}
	}
	
}