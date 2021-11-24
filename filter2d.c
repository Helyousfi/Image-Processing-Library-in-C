
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <filter2d.h>


void filter2d_mean3(double * img_in, long width, long height, double * img_out)
{
	double Hi = 1.0 / 9.0;
	double H[3][3] = {{ Hi, Hi, Hi},
					  { Hi, Hi, Hi},
					  { Hi, Hi, Hi}};
	int M = 1;
	int N = 1;
	int x, y;

	double s;


	for (int i = 0; i < height * width; i++)
		img_out[i] = img_in[i];

	for (y = 1; y < height-1; y++)
	{
		for (x = 1; x < width-1; x++)
		{
			s = 0.0;
			for (int m = -M; m <= M; m++)
			{
				for (int n = -N; n <= N; n++)
				{
					s += img_in[(y - m)*width + x - n] * H[m+1][n+1];			
				}
			}
			img_out[y*width + x] = s;
		}
	}
}


void sobel_norm(double * img_in, long width, long height, double * img_out)
{
	double Sx[3][3] = { { 1.0/8.0, 0, -1.0/8.0},
					  { 2.0/8.0, 0, -2.0/8.0},
					  { 1.0/8.0, 0, -1.0/8.0} };
	double Sy[3][3] = { { 1.0/8.0, 2.0/8.0, 1.0/8.0},
					  { 0.0, 0.0, 0.0},
					  { -1.0/8.0, -2.0/8.0, -1.0/8.0} };


	int M = 1;
	int N = 1;
	int x, y;

	double sx, sy;


	for (int i = 0; i < height * width; i++)
		img_out[i] = 0;

	for (y = 1; y < height - 1; y++)
	{
		for (x = 1; x < width - 1; x++)
		{
			sx = sy = 0.0;
			for (int m = -M; m <= M; m++)
			{
				for (int n = -N; n <= N; n++)
				{
					sx += img_in[(y - m) * width + x - n] * Sx[m + 1][n + 1];
					sy += img_in[(y - m) * width + x - n] * Sy[m + 1][n + 1];
				}
			}
			img_out[y * width + x] = sqrt(sx*sx+sy*sy);
		}
	}
}



static double * gaussian2d_create(double sigma, long * size)
{
	// calcul de taille
	int t;
	double * g = NULL, s = 0.0, v;

	t = 1+2* (int)ceil(3 * sigma);
	*size = t;

	// allocation
	g = (double*)calloc(t * t, sizeof(double));

	
	
	//remplissage et calcul de somme
	for(int j = 0; j < t; j++)   //La deuxième lettre c la verticale
	{
		double y = (double)(j - t/2) / sigma;
		for (int i = 0; i < t; i++)
		{
			double x = (double)(i - t/2) / sigma;
			v = exp(-(x * x + y * y) / 2.0);
			s += v;
			g[t * j + i] = v;
		}
	}


	// normalisation


	/*
	double sum = 0.0;
	for (int i = 0; i < t * t; i++)
	{
		sum = sum + g[i];
	}
	*/
	

	for (int k = 0; k < t * t; k++)
	{
		g[k] /= s;
	}
	

	return g;
}




static void convolution2d(double * img_in, long width, long height,
									        double * mask2d, long tx, long ty,
									        double * img_out)
{
	for (int i = 0; i < height * width; i++)
		img_out[i] = img_in[i];

	double s;
	for (int y = ty/2; y < height - ty/2; y++)
	{
		for (int x = tx/2; x < width - tx/2; x++)
		{
			s = 0.0;
			for (int m = -ty/2; m <= ty/2; m++)
			{
				for (int n = -tx/2; n <= tx/2; n++)
				{
					s += img_in[(y - m) * width + x - n] * mask2d[(m + ty / 2) * tx + n + tx / 2];
				}
			}
			img_out[y * width + x] = s;
		}
	}


	// for (int m = 0; m < ty; m++)
	// s += img_in[(y - (m-ty/2)) * width + x - (n+tx/2)] * mask2d[m * tx + n];
	// s += img_in[(y + (m-ty/2)) * width + x + (n+tx/2)] * mask2d[(ty-1-m) * tx + tx - 1 - n];
}






static double ** gaussian2d_create_matrix(double sigma, long * size)
{
	int t, j;
	double** g = NULL, s = 0.0, v;

	t = 1 + 2 * (int)ceil(3 * sigma);
	*size = t;

	// allocation
	g = (double**)calloc(t, sizeof(double *));
	for (j = 0; j < t; j++) g[j] = (double*)calloc(t, sizeof(double));
	
	//!
	g[0] = (double*)calloc(t * t, sizeof(double));
	for (j = 1; j < t; j++) g[i] = g[0] + j * t;
	//! 

	//remplissage et calcul de somme
	for (int j = 0; j < t; j++)  
	{
		double y = (double)(j - t / 2) / sigma;
		for (int i = 0; i < t; i++)
		{
			double x = (double)(i - t / 2) / sigma;
			v = exp(-(x * x + y * y) / 2.0);
			s += v;
			g[j][i] = v;
		}
	}


	for (int k = 0; k < t * t; k++)
	{
		g[k] /= s;
	}

	return NULL;
}




static void convolution2d_by_matrix(double * img_in, long width, long height,
													  double ** mask2d, long tx, long ty,
									                  double * img_out)
{
	for (int i = 0; i < height * width; i++)
		img_out[i] = 0.0;

	double s;
	for (int y = ty / 2; y < height - ty / 2; y++)
	{
		for (int x = tx / 2; x < width - tx / 2; x++)
		{
			s = 0.0;
			for (int m = -ty / 2; m <= ty / 2; m++)
			{
				for (int n = -tx / 2; n <= tx / 2; n++)
				{
					s += img_in[(y - m) * width + x - n] * mask2d[l][k];
				}
			}
			img_out[y * width + x] = s;
		}
	}
}


void filter2d_gaussian(double * img_in, long width, long height,
											 double sigma,
											 double * img_out)
{
	double* mask2d = NULL;
	int t; 
	// = (int*)calloc(1, sizeof(int));
	
	mask2d = gaussian2d_create(sigma, &t);

	convolution2d(img_in, width, height, mask2d, t, t, img_out);
	free(mask2d);

	double** filter = NULL;
	int t, j;

	filter = gaussian2d_create_matrix(sigma, &t);
	convolution2d_by_matrix(img_in, width, height, filter, t, t, img_out);
	for (j = 0; j < t; j++) free(filter[j]);
	free(filter);

	free(filter[0]);
	free(filter);
}


void img_get_raw(double * img, long width, long height,
		             long no,
						     double * v)
{
	for (int i = 0; i < width; i++)
		v[i] = img[no * width + i];
}


void img_set_raw(double * img, long width, long height,
								 long no,
								 double * v)
{
	for (int i = 0; i < width; i++)
		img[no * width + i] = v[i];
}


void img_get_column(double * img, long width, long height,
									  long no,
										double * v)
{
	for (int j = 0; j < height; j++)
		v[j] = img[j * width + no];
}


void img_set_column(double * img, long width, long height,
										long no,
										double * v)
{
	for (int j = 0; j < height; j++)
		img[j * width + no] = v[j];
}


static double * gaussian1d_create(double sigma, long * size)
{
	// calcul de taille
	int t;
	double* g = NULL, s = 0.0, v;

	t = 1 + 2 * (int)ceil(3.0 * sigma);
	*size = t;

	// allocation
	g = (double*)calloc(t, sizeof(double));


	//remplissage et calcul de somme

	for (int i = 0; i < t; i++)
	{
		double x = (double)(i - t / 2) / sigma;
		v = exp(-(x * x) / 2.0);
		s += v;
		g[i] = v;
	}
	


	// normalisation

	for (int k = 0; k < t; k++)
	{
		g[k] /= s;
	}


	return g;
}


static void convolution1d(double * v_in, long size,
													double * mask1d, long t,
													double * v_out)
{
	for (int i = 0; i < size; i++)
		v_out[i] = 0.0;

	double s;
	
	for (int x = t / 2; x < size - t / 2; x++)
	{
		s = 0.0;
		for (int k=0; k < t; k++)
		{
			s += v_in[x-(k-t/2)] * mask1d[k];
		}
		v_out[x] = s;
	}
	
}


void filter2d_gaussian_fast(double * img_in, long width, long height,
														double sigma,
														double * img_out)
{
	double* v_in = NULL, * v_out = NULL, * filter = NULL;
	int i, j, t;

	filter = gaussian1d_create(sigma, &t);

	v_in = (double*)calloc(width, sizeof(double));
	v_out = (double*)calloc(width, sizeof(double));

	for (j = 0; j < height; j++)
	{
		
		img_get_raw(img_in, width, height, j, v_in);
		convolution1d(v_in, width, filter, t, v_out);
		img_set_raw(img_out, width, height, j, v_out);
	}
	

	for (i = 0; i < height; i++)
	{

		img_get_column(img_out, width, height, i, v_in);
		convolution1d(v_in, width, filter, t, v_out);
		img_set_column(img_out, width, height, i, v_out);
	}

	free(v_in);
	free(v_out);
	free(filter);
}

static double v_min(double* v, long size)
{
	double temp_min = v[0];
	for (int i = 0; i < size; i++)
	{
		if (v[i] < temp_min)
		{
			temp_min = v[i];
		}
	}
	return temp_min;
}
static double v_max(double* v, long size)
{
	double temp_max = v[0];
	for (int i = 0; i < size; i++)
	{
		if (v[i] > temp_max)
		{
			temp_max = v[i];
		}
	}
	return temp_max;
}

static double v_mean(double* v, long size)
{
	double temp_mean = 0.0;
	for (int i = 0; i < size; i++)
	{
		temp_mean += v[i];
	}
	return temp_mean/(double)size;
}



int compare(const void* v1, const void* v2)
{
	double* e1 = ((double*)v1);
	double* e2 = ((double*)v2);

	// 0 = égalité, >0 "v1" > "v2", <0 "v1" < "v2"
	if (*e1 > *e2)
		return 1;
	return *e1 < *e2;
}

static double v_median(double* v, long size)
{
	qsort((void*)v, (size_t)size, sizeof(double), compare);
	return v[size/2];
}



typedef double (* PROC)(double *, long);


static void filter2d_generic(double * img_in, long width, long height,
							  						 long tx, long ty, PROC proc,
														 double * img_out)
{
	int i, x, y, k, l;
	double s, * v = NULL;
	for (int i = 0; i < height * width; i++) img_out[i] = img_in[i];
	v = (double*)calloc(tx * ty, sizeof(double));
	int a = 0;
	for (y = ty / 2; y < height - ty / 2; y++)
		for (x = tx / 2; x < width - tx / 2; x++)
		{
			//Remplissage de v
			i = 0;
			for (l = -ty/2; l < ty/2; l++)
			{
				for (k = -tx / 2; k <= tx / 2; k++)
				{
					v[i++] = img_in[(y + l) * width + x + k];
				}
			}
			img_out[y * width + x] = proc(v, tx*ty);
		}
	free(v);
}


void filter2d_by_method(double * img_in, long width, long height,
	                      long tx, long ty, long method,
												double * img_out)
{
	PROC v_proc[4] = { v_min, v_max, v_mean, v_median };
	filter2d_generic(img_in, width, height, tx, ty, v_proc[method], img_out);
}

typedef double (*PROC1)(double*, long);
typedef int (*PROC2)(int*, long);
typedef double (*PROC_NULL)(void);

proc3 = (PROC_NULL *)proc1;
proc4 = (PROC1 *)proc3;




