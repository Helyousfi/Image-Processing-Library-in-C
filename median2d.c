#include <stdlib.h>

#include <median2d.h>



int compare(const void * v1, const void * v2)
{
	double* e1 = ((double*)v1);
	double* e2 = ((double*)v2);

	// 0 = égalité, >0 "v1" > "v2", <0 "v1" < "v2"
	if (*e1 > *e2)
		return 1;
	return *e1 < *e2;

	/*
	double e1 = *((double*)v1);
	double e2 = *((double*)v2);

	// 0 = égalité, >0 "v1" > "v2", <0 "v1" < "v2"
	if (e1 > e2)
		return 1;
	else if (e1 < e2)
		return -1;
	return 0;
	*/

	/*
	if (*((double*)v1) < *((double*)v2))
		return -1;
	return  (*((double*)v1) > *((double*)v2));
	*/
	
}




void filter2d_median(double * img_in, long width, long height,
										 long tx, long ty,
										 double * img_out)
{
	int i, x, y, k, l;
	double s, *v=NULL;
	for (int i = 0; i < height * width; i++) img_out[i] = img_in[i];
	
	int a = 0;
	for (y = ty / 2; y < height - ty / 2; y++)
		for (x = tx / 2; x < width - tx / 2; x++)
		{
			//Remplissage de v
			i = 0;
			for (l = -ty / 2; l < ty / 2; l++) // for (l = -ty / 2; l <= ty / 2; l++)
			{
				for (k = -tx / 2; k <= tx / 2; k++)
				{
					v[i++] = img_in[(y + l) * width + x + k];
				}
			}

			//tri
			qsort((void*)v, (size_t)(tx * ty), sizeof(double), compare);
			
			// sortie = valeur "milieu" du tri
			
			img_out[y * width + x] = v[tx * ty / 2];
		}


	free(v);
}