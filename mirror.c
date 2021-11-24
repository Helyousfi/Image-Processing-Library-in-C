
#include <mirror.h>


void img_mirror_horizontal(double * img_in, long width, long height,
													 double * img_out)
{

	
	/*
	for (int j = 0; j < height; j++)
		for (int i = 0; i < width; i++)
			img_out[width-1 - i + width * j] = img_in[i + width * j];


	for (int j = 0; j < height; j++)
		for (int i = 0; i < width; i++)
			*(img_out-- + (width - 1) * (j + 1)) = *(img_in++);

	*/
	

	img_out += width - 1;
	for (int j = 0; j < height; j++)
		for (int i = 0; i < width; i++)
			*(img_out--) = *(img_in++);
		img_out += 2*width;
	

}