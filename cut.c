

#include <cut.h>
#include <cstddef>

typedef struct _RGBA_
{
	unsigned char r;
	unsigned char g;
	unsigned char b;
	unsigned char a;
} RGBA;

void img_extract_area(unsigned char * img_color_in, long width, long height,
											long x, long y, long dx, long dy,
											unsigned char * img_color_out)
{
	/*
	// Version "classique" unsigned char
	for (int j = 0; j < dy; j++)
		for (int i = 0; i < dx; i++)
			for (int c = 0; c < 4; c++)
				img_color_out[(i + dx * j) * 4 + c] = img_color_in[(i + x + c + (j + y) * width) * 4];

	
	
	// Version "classique" struct RGBA
	RGBA *img_rgba_in = (RGBA *)img_color_in, * img_rgba_out = (RGBA*)img_color_out;
	for (int j = 0; j < dy; j++)
		for (int i = 0; i < dx; i++)
			img_rgba_out[j * dx + i] = img_rgba_in[(j + y) * width + i + x];
	*/
	// Arithmetique des pointeurs
	RGBA* img_rgba_in = (RGBA*)img_color_in, * img_rgba_out = (RGBA*)img_color_out;

	img_rgba_in += y * width + x;
	for (int j = 0; j < dy; j++)
	{
		for (int i = 0; i < dx; i++)
			*(img_rgba_out++) = *(img_rgba_in++);
	}
	img_rgba_in += width - dx;
}


#define MULADD(a, b, c) (((a)*(b))+(c)) // PARENTHESIS ARE IMPORTANT!!
#define MAX(a,b) (a > b ? a : b)


static int compare(const void* v1, const void* v2)
{
	RGBA* e1 = (RGBA*)v1;
	RGBA* e2 = (RGBA*)v2;

	unsigned char m1, m2;

	m1 = MAX(e1->r, MAX(e1->g, e1->b));
	m2 = MAX(e2->r, MAX(e2->g, e2->b));
	if (m1 > m2) return 1;
	if (m1 < m2) return -1;

	/*
	if (e1->r > e1->g)
	{
		if (e1->r > e1->b)
		{
			m1 = e1->r;
		}
		else
		{
			m1 = e1->b;
		}
	}
	else
	{
		if (e1->g > e1->b)
		{
			m1 = e1->g;
		}
		else
		{
			m1 = e1->b;
		}

	}
	*/
	return 0;
}



void img_diffuse_hot_color(unsigned char* img_color_in,
							long width, long height,
							long tx, long ty,
							unsigned char* img_color_out)
{
	int i, x, y, k, l;
	RGBA * img_rgba_in = (RGBA *)img_color_in, 
		 * img_rgba_out = (RGBA *)img_color_out,
		 * v = NULL;

	for (int i = 0; i < height * width; i++) img_color_out[i] = img_color_in[i];

	int a = 0;
	v = (RGBA*)calloc(tx * ty, sizeof(RGBA));
	for (y = ty / 2; y < height - ty / 2; y++)
		for (x = tx / 2; x < width - tx / 2; x++)
		{
			//Remplissage de v
			i = 0;
			for (l = ty/2; l < ty; l++)
			{
				for (k = -tx / 2; k <= tx / 2; k++)
				{
					v[i++] = img_rgba_in[(y + l) * width + x + k];
				}
			}

			//tri
			qsort((void*)v, (size_t)(tx * ty), sizeof(RGBA), compare);

			// sortie = valeur "milieu" du tri

			img_rgba_out[y * width + x] = v[tx * ty / 2];
		}


	free(v);
}