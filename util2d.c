#include <stdio.h>
#include "util2d.h"
#include <math.h>

void print_matrix(const char *name, double **inmat, int n1, int n2)
{
     int i,j;
     printf("\n contents of %s : \n", name);
     for( i = 0; i < n1; i++)
     {
	  for( j = 0; j < n2; j++)
	       printf("%1.16f, ", inmat[i][j]);
	  printf("\n"); 

     }
     //done!
}

void print_array(const char *name, double *inmat, int n1)
{
     int i;
     printf("\n contents of %s : \n", name);
     for( i = 0; i < n1; i++)
	  printf("%1.16f \n", inmat[i]);
     //done!
}

double max_abs_array(double *inmat, int n1)
{
     int i;
     double max = fabs(inmat[0]);
     for( i = 1; i < n1; i++)
	  if(fabs(inmat[i]) > max)
	       max = fabs(inmat[i]);

     return max;
     //done!
}
