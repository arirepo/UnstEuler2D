#include <stdio.h>
#include <stdlib.h>
#include "solve_lu_serial.h"

/*  / ____|/ __ \| |\ \    / /  ____|/ ____| */
/* | (___ | |  | | | \ \  / /| |__  | (___   */
/*  \___ \| |  | | |  \ \/ / |  __|  \___ \  */
/*  ____) | |__| | |___\  /  | |____ ____) | */
/* |_____/ \____/|______\/   |______|_____/  */
/* /\  _  \                         /\ \ */
/* \ \ \ \ \   __  _                \ \ \____ */
/*  \ \  __ \ /\ \/ \      _______   \ \  __ \ */
/*   \ \ \/\ \\/>  </     /\______\   \ \ \ \ \ */
/*    \ \_\ \_\/\_/\_\    \/______/    \ \_ __/ */
/*     \/_/\/_/\//\/_/                  \/___/ */
/*           (_)                               */
/*  _   _ ___ _ _ __   __ _                    */
/* | | | / __| | '_ \ / _` |                   */
/* | |_| \__ \ | | | | (_| |                   */
/*  \__,_|___/_|_| |_|\__, |                   */
/*                     __/ |                   */
/*                    |___/                    */
/*                    /__/\                    */
/*                    \  \:\                   */
/*   ___     ___       \  \:\                  */
/*  /__/\   /  /\  ___  \  \:\                 */
/*  \  \:\ /  /:/ /__/\  \__\:\                */
/*   \  \:\  /:/  \  \:\ /  /:/                */
/*    \  \:\/:/    \  \:\  /:/                 */
/*     \  \::/      \  \:\/:/                  */
/*      \__\/        \  \::/                   */
/*                    \__\/                    */

/* 
   task : solves Ax=b using LU decomposition already stored in A 
   and the permutation matrix P. Use lu_serial.c to do LU decomposition first. 
   The size of A is nxn. 
   ver : 1.00
   ghasemi.arash@gmail.com
*/
int solve_lu_serial(int *P, double *A, double *x, double *b, int n)
{
  //allocating local variables 
  int i,j;
  double *Ap = (double *)malloc(n*n*sizeof(double));
  double *bp = (double *)malloc(n*sizeof(double));
  double *x_bar = (double *)malloc(n*sizeof(double));

  //transforming A,b to fully pivotted Ap = P*A=[L][U] and bp = P*b
  // so final solution [x] is [Ap][x]=[bp] -> [L][U] [x] = [bp]
  // step 1 : forward substitution [U] [x] = [L]^-1 [bp] = [x_bar]
  // step 2: backward substitution [x] = [U]^-1 [x_bar]
  
  serial_matrix_mult(P, A, Ap, n, n , n);
  serial_matrix_mult(P, b, bp, n, n , 1);
  //starting forward substitution
  for(i = 0; i < n; i++)
    {
      x_bar[i] = bp[i];
      for(j = 0; j < i; j++)
	x_bar[i] -= Ap[i*n+j] * x_bar[j];
    }  
  //starting backward substitution
  for(i = (n-1); i >= 0; i--)
    {
      x[i] = x_bar[i]/Ap[i*n+i];
      for(j = (n-1); j > i; j--)
	x[i] -= (Ap[i*n+j]/Ap[i*n+i]) * x[j];
    }  
  //clean up
  free(Ap);
  free(bp);
  free(x_bar);

  //completed successfully!
  return 0; 
}
/* the following performs a serial matrix multipication C = A*B
   where A=A(m,n) and B=B(n,l) hence C=C(m,l).
   NOTE: IT DOESN'T CHECK FOR DIMENSIONAL COMPATIBILITY
   SO SPECIAL CARE MUST BE TAKEN BEFORE USING THIS.
*/
void serial_matrix_mult(int *A, double *B,double *C,int m, int n , int l)
{
     /* local variables */
     int i,j,k;
     for ( i= 0; i < m ; i++)
	  for ( j= 0; j < l ; j++)
	  {
	       C[i*l + j] = 0.;
	       for ( k= 0; k < n ; k++)	  
		 C[i*l + j] += ((double)A[i*n+k]) * B[k*l+j];
	  }
     /* simply finished .*/     
}
