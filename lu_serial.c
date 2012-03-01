#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lu_serial.h"

	/* +------------------------------------------------------------+ */
        /* |..____  _____ ____..___...._...._..............._...._..._..| */
	/* |./ ___|| ____|  _ \|_ _|../ \..| |.............| |..| |.| |.| */
	/* |.\___ \|  _| | |_) || |../ _ \.| |......_____..| |..| |.| |.| */
	/* |..___) | |___|  _ < | |./ ___ \| |___..|_____|.| |__| |_| |.| */
	/* |.|____/|_____|_|.\_\___/_/...\_\_____|.........|_____\___/..| */
	/* |............................................................| */
	/* |.........Arash.Ghasemi.-.ghasemi.arash@gmail.com............| */
	/* +------------------------------------------------------------+ */
/* 
   task : performs serial LU-decomposition of [A] and stores the results in [A] 
   itself. The permuation matrix [P] designates the correct ordering of rows.
   ver : 1.00 
*/

int lu_serial(double *A, int *P, int n)
{
  //local vars in stack
  int i,j,k;
  int max_i = 0;

  // allocating and initializing the permutation array p
  int *p = (int *)malloc(n * sizeof(int));
  for( i = 0; i < n; i++)
    p[i] = i;

  //main loop - finishes when LU is complete!
  for( k = 0 ; k < (n-1); k++)
    {
      //step 1 - partial pivoting 
      // first find max value in the column
      find_max(A, p, k, n, &max_i); 
      if( p[max_i] != k) //then pivot the row
	pivot_p(p, k, max_i); //swaps k, max_i
      //compute multipliers column m
      for (i = (k+1); i < n; i++)
	A[p[i]*n+k] = A[p[i]*n+k] / A[p[k]*n+k];
      //updating As
      for (i = (k+1); i < n; i++)
	for(j = (k+1); j < n; j++)
	  A[p[i]*n+j] -= (A[p[i]*n+k] * A[p[k]*n+j]);
    }

  //setting permutation matrix
  // first resetting the input matrix
  for( i = 0; i < n; i++)
    for( j = 0; j < n; j++)
      P[i*n+j] = 0;
  //then setting it up!
  for( i = 0; i < n; i++)
      P[i*n+p[i]] = 1;
 
  //little clean-up
  free(p);

  //completed successfully!
  return 0;

}

inline int find_max(double *A, int *p, int k, int n, int *max_i)
{
  int i;
  double max_val = 0., this_val = 0.; 
  for (i = k; i < n; i++)
    {
      this_val = fabs(A[p[i]*n+k]);
      if( this_val > max_val)
	{
	  max_val = this_val;
	  *max_i = i;
	}
    }
 
  //completed successfully!
  return 0;

}

inline int pivot_p(int *p, int i, int j)
{
  int temp = p[j];
  p[j] = p[i];
  p[i] = temp; 

  //completed successfully!
  return 0;
}  

void print_matrix_double(const char *name, double *inmat, int n1, int n2)
{
     int i,j;
     printf("\n contents of %s : \n", name);
     for( i = 0; i < n1; i++)
     {
	  for( j = 0; j < n2; j++)
	       printf("%e, ", inmat[i*n2+j]);
	  printf("\n"); 

     }
     //done!
}

void print_matrix_int(const char *name, int *inmat, int n1, int n2)
{
     int i,j;
     printf("\n contents of %s : \n", name);
     for( i = 0; i < n1; i++)
     {
	  for( j = 0; j < n2; j++)
	       printf("%d, ", inmat[i*n2+j]);
	  printf("\n"); 

     }
     //done!
}

void print_array_double(const char *name, double *inmat, int n1)
{
     int i;
     printf("\n contents of %s : \n", name);
     for( i = 0; i < n1; i++)
	  printf("%e \n", inmat[i]);
     //done!
}
