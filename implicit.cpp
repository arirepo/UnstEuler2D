#include <stdio.h>
#include "implicit.h"
#include "maps.h"
#include <vector>
#include <algorithm>
#include "residuals.h"
#include "util2d.h"
#include "gauss_seidel_valid.h"

using namespace std;

//allocates sparse matrices A and b required for left and right hand sides of implicit temporal discretization.
//      INPUTS
// nn : number of nodes
// neqs : number of equations
// :
//      OUTPUTS
// nnz : total number of non-zero blocks in the sparse matrix
// *ia : array containing the zero based start of each raw
// *ja : array containing the zero based column number for row[ia[i] - ia[i+1]-1] 
// *iau : array containing the zero based position of diagonal for row i
// *A : one-dimensional array containing all two-dimensional sub blocks in sparse form
// *rhs: one-dimensional array containing rhs 
int alloc_A_b( int nn, int neqs, int nt, int **tri_conn, int *nnz, int **ia, int **ja,  int **iau, double **A, double **rhs)
{
  //local indices
  int i, j, k; 
  vector<int>::iterator it; //an iterator

  //allocating p2p map ...
  vector<int> *p2p = (vector<int> *)calloc( nn , sizeof(vector<int>));
  // creating p2p map ...
  create_p2p(nt, tri_conn, nn, p2p);
  //debugging
  // for( i = 0; i < nn; i++)
  //   { 
  //     printf("\n p2p[%d]: ", i);
  //     fflush(stdout);
  //     for ( j = 0; j < p2p[i].size(); j++)
  // 	printf("%d,", p2p[i][j]);
  //   }
  // printf("\n");

  //calculating nnz
  *nnz = 0; //hard reset!
  for ( i = 0 ; i < nn; i++)
    (*nnz)+= p2p[i].size();

  // printf("\n nnz: %d", *nnz);

  //allocating sparse-matrix arrays using our p2p map
  //THESE WILL BE RETURNED BACK
  (*ia) = (int *)malloc( (nn + 1) * sizeof(int) );
  (*iau) = (int *)malloc( nn * sizeof(int) );
  (*ja) = (int *)malloc( (*nnz) * sizeof(int) );
  (*A) = (double *)calloc( neqs * neqs * (*nnz) , sizeof(double) );
  (*rhs) = (double *)calloc( neqs * nn , sizeof(double) );

  // filling ja
  k = 0;
  for ( i = 0; i < nn; i++)
    for( j = 0; j < p2p[i].size(); j++) 
      {
	(*ja)[k++] = p2p[i][j];
	// printf("\n (*ja)[k]: %d", (*ja)[k-1]);
      }

  // filling ia, aiu
  (*ia)[0] = 0;
  (*iau)[0] = 0;
  for ( i = 1 ; i < nn; i++)
    {
      (*ia)[i] = (*ia)[i-1] + p2p[i-1].size();
      // printf("\n (*ia)[i]: %d", (*ia)[i]);
      //find diagonal element
      it = find (p2p[i].begin(), p2p[i].end(), i);
      (*iau)[i] = (*ia)[i] + (int)(it - p2p[i].begin());      
      // printf("\n (*iau)[i]: %d", (*iau)[i]);
    }
  (*ia)[i] = (*ia)[i-1] + p2p[i-1].size()-1;
  // printf("\n (*ia)[i]: %d", (*ia)[i]);

  // clean - up
  free(p2p);

  //completed successfully!
  return 0;

}

//implemets euler explicit scheme in Ax = b = rhs form
int Axb_euler_explicit(double *Q, double *Q_inf, double gamma, double CFL, int ITR_MAX, int itr_per_msg, double *x, double *y, int *bn_nodes, int nn, int neqs, int nt, int **tri_conn, int nnz, int *ia, int *ja,  int *iau, double *A, double *rhs)
{

  double *xn = (double *)malloc( neqs * nn * sizeof(double) );
  double *xn1 = (double *)malloc( neqs * nn * sizeof(double) );
  double *x_star = (double *)malloc( neqs * nn * sizeof(double) );
  //locals
  int i,j,k;
  double *R = (double *)calloc((nn*neqs) , sizeof(double));
  int ITR = 0;
  double *int_uplusc_dl = (double *)calloc(nn , sizeof(double) );

  // main iteration loop
  for( ITR = 1; ITR <= ITR_MAX; ITR++)
    {

      //finding the residuals
      calc_residuals( Q, Q_inf, gamma, nn, neqs, x, y, nt, tri_conn, bn_nodes, R);

      // calculating line integral int( (|u_bar| + c) dl ) 
      calc_int_uplusc_dl( Q, gamma, neqs, nn, x, y, nt, tri_conn, bn_nodes, int_uplusc_dl);

      for ( i = 0; i < nn; i++)
	for( j = 0; j < neqs; j++)
	  rhs[i*neqs + j] = -R[i*neqs + j];

      for ( i = 0; i < nn; i++)
	{
	  k = iau[i];
	  for( j = 0; j < neqs; j++)
	    A[k*neqs*neqs + j*neqs+j] = int_uplusc_dl[i]/CFL;
	}

      gauss_seidel_solve_pivoting(nn, nnz, ia, ja,  iau, A, rhs, neqs, x_star, xn1, xn);

      //updating Q
      for( i = 0; i < nn; i++)
	for( j = 0; j < neqs; j++)
	  Q[i*neqs + j] = Q[i*neqs + j] + xn[i*neqs + j];


      if(!(ITR % itr_per_msg)) // show status each itr_per_msg time
	printf("ITR = %d, max_abs(R[1,2,3,4]) = %17.17e\n", ITR, max_abs_array(xn, (neqs*nn)));


    }

  printf("The infinity condition is:\n");
  print_array("Q_inf", Q_inf, 4);

  //clean-up
  free(R);
  free(int_uplusc_dl);
  free(xn);
  free(xn1);
  free(x_star);



  //completed successfully!
  return 0;

}
