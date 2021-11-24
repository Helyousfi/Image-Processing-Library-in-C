

#include <stretch.h>


static void v_compute_min_max(double * v_in, long size, double * v_min, double * v_max)
{
	int i;
	double v, vmin, vmax;

	vmin = vmax = v_in[0];

	for (i = 1; i < size; i++)
	{
		v = v_in[i];
		if (v < vmin) {
			vmin = v;
		}
		else if (v > vmax) {
			vmax = v;
		}
	}

	*v_min = vmin;
	*v_max = vmax;
}


void img_stretch_intensity(double * img_in, long width, long height, double * img_out)
{
	double v_min, v_max, a1 = 0.0, a2 = 255.0;
	int i, size = width * height;

	v_compute_min_max(img_in, size, &v_min, &v_max);
	if (v_max != v_min)
	{
		for (i = 0; i < size; i++)
		{
			img_out[i] = a1 + (a2 - a1) / (v_max - v_min) * (img_in[i] - v_min);
		}
	}
	else
	{
		for (i = 0; i < size; i++)
		{
			img_out[i] = (a2 + a1) / 2.0;
		}
	}
		
	


}