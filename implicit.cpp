#include <stdio.h>
#include "implicit.h"
#include "maps.h"
#include <vector>
#include <algorithm>
#include "residuals.h"
#include "util2d.h"
#include "gauss_seidel_valid.h"
#include "flux.h"

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

//accumulate Df sub-block to the diagonal of matrix A located at (row, row)
inline int accumulate_to_diag(double *A, double **Df, int row, int neqs, int *iau)
{
  //locals
  int i,j;
  for( i = 0; i < neqs; i++)
    for( j = 0; j < neqs; j++)
      A[iau[row]*neqs*neqs + i * neqs + j] += Df[i][j];

  //completed successfully!
  return 0;
}

inline void identity(double *S, int neqs)
{
  int i,j;
  for( i = 0; i < neqs; i++)
    for( j = 0; j < neqs; j++)
      if ( i == j)
	S[i*neqs + j] = 1.;
      else
	S[i*neqs + j] = 0.;

}
inline void freez(double *S, int n1, int n2 )
{
  int i,j;
  for( i = 0; i < n1; i++)
    for( j = 0; j < n2; j++)
      S[i*n2 + j] = 0.;


}		    

//accumulate Df sub-block to the locaation (i,j) of matrix A.
inline int accumulate_to_ij(double *A, double **Df, int i, int j, int neqs, int *ia, int *ja)
{

  int l, m, loc;
  int i_start = ia[i];
  int i_end = ia[i+1]-1;
  loc = i_start; //initial location
  for( l = i_start; l <= i_end; l++) //adding appropriate offset
    if(ja[l] == j)
      break;
    else
      loc++;  
  // now we have the location of the block which we should accumulate to.
  for( l = 0; l < neqs; l++)
    for( m = 0; m < neqs; m++)
      A[loc*neqs*neqs + l * neqs + m] += Df[l][m];


  //completed successfully!
  return 0;
}


// first reset A and b
// fills  martices A and b with Jacobians and residuals respectively
int fill_A_b(double *Q, double *Q_inf, double gamma, double *x, double *y, int *bn_nodes, int nn, int neqs, int nt, int **tri_conn, int nnz, int *ia, int *ja,  int *iau, double *A, double *rhs)
{

     //local vars
  int i, j, t, kk;
     int n_right, n_left;
     double xc, yc, xmid, ymid;
     double nx, ny;

     double *n_hat = (double *)calloc(2 , sizeof(double));
     int *f_select = (int *)calloc(4 , sizeof(int));
     
     double *fvl_p = (double *)calloc( neqs , sizeof(double));
     double *fvl_m = (double *)calloc( neqs , sizeof(double));
     double *Q_edge = (double *)calloc( neqs , sizeof(double));
     double *fw = (double *)calloc( neqs , sizeof(double)); //wighted wall flux
  
     //allocating Van Leer flux vector jacobians d_+- 
     double **d_fvl_p = (double **)calloc( neqs , sizeof(double *));
     double **d_fvl_m = (double **)calloc( neqs , sizeof(double *));
     double **dfw = (double **)calloc( neqs , sizeof(double *));//wighted wall flux Jaco

     for( i = 0; i < neqs; i++)
     {
	  d_fvl_p[i] = (double *)calloc( neqs , sizeof(double));
	  d_fvl_m[i] = (double *)calloc( neqs , sizeof(double));
	  dfw[i] = (double *)calloc( neqs , sizeof(double));
     }

     //resseting the A
     for( i = 0; i < nnz; i++)	  
       for(j=0; j < neqs; j++)
	 for(t=0; t < neqs; t++)
	   A[neqs*neqs*i+neqs*j+t] = 0.;

     //resseting the rhs vector
     for( i = 0; i < nn; i++)	  
	  for(j=0; j < neqs; j++)
	       rhs[neqs*i+j] = 0.;
     
     //main loop : loop over all triangles
     for (t = 0; t < nt; t++)
     {
	  //calculate the center of triangle - mean average- 
	  xc = (x[tri_conn[t][0]] + x[tri_conn[t][1]] + x[tri_conn[t][2]]) / 3.;
	  yc = (y[tri_conn[t][0]] + y[tri_conn[t][1]] + y[tri_conn[t][2]]) / 3.;
	  
	  for(i = 0; i < 3; i++) //loop over vertices of each triangle
	  {
	       //determine left and right
	       n_left = i;	       
	       n_right = (i<2)?(i+1):0;

	       //converting local to global node number 
	       n_left = tri_conn[t][n_left];
	       n_right = tri_conn[t][n_right];

	       if((bn_nodes[n_left] == 1) && (bn_nodes[n_right] == 1)) //both nodes are on the free stream boundaries.
	       {
		 //right now skip, because we fix DQ = 0.
		 //we will come back and improve it for sure.   
		    
	       }
	       if((bn_nodes[n_left] > 1) && (bn_nodes[n_right] > 1)) //both nodes are on the wall boundaries.
	       {
		    nx = 0.5*(y[n_right] - y[n_left]);
		    ny = -0.5*(x[n_right]-x[n_left]);
		    n_hat[0] = nx;
		    n_hat[1] = ny;

		    // IMPORTANT: To store values of wall fluxes which are different from van-leer fluxes, we still use same buffer buffers fvl_m and d_fvl_m which are originally intended to be used for storing left and right van-leer fluxes f+ and f-. This reduces the extra programming and memory usage needed for additional arrays with different names.
		    ///    (fvl_m,d_fvl_m)o-----wall--------o (fvl_p,d_fvl_p)

		    //contributing to the left node 
		    calc_wall_flux((Q+n_left*neqs), fvl_m, d_fvl_m, neqs, gamma, n_hat);

		    //contributing to the right node 
		    calc_wall_flux((Q+n_right*neqs),fvl_p, d_fvl_p, neqs, gamma, n_hat);

		    //weighted averages of flux and jacobians for the left node
		    for(j=0; j < neqs; j++)
			 fw[j] = .75*fvl_m[j] + .25*fvl_p[j];

		    for(j=0; j < neqs; j++)
		      for(kk=0; kk < neqs; kk++)
			dfw[j][kk] = .75*d_fvl_m[j][kk] + .25*d_fvl_p[j][kk];

		    //accumulate wall flux with negative sign to the residual of the left node
		    for(j=0; j < neqs; j++)
		    {
			 rhs[neqs*n_left+j] -= fw[j];
		    }

		    //accumulate Jacobian to the diagonal of matrix A for left nodes
		    accumulate_to_diag(A,dfw,n_left,neqs, iau);

		    //weighted averages of flux and jacobians for the right node
		    for(j=0; j < neqs; j++)
			 fw[j] = .25*fvl_m[j] + .75*fvl_p[j];

		    for(j=0; j < neqs; j++)
		      for(kk=0; kk < neqs; kk++)
			dfw[j][kk] = .25*d_fvl_m[j][kk] + .75*d_fvl_p[j][kk];

		    //accumulate wall flux with negative sign to the residual of the right node
		    for(j=0; j < neqs; j++)
		    {
			 rhs[neqs*n_right+j] -= fw[j];
		    }

		    //accumulate Jacobian to the (n_left, n_right) of matrix A 
		    accumulate_to_ij(A, dfw, n_left, n_right, neqs, ia, ja);
		    
	       }

	       //calculating the center of that edge
	       xmid = (x[n_right] + x[n_left]) / 2.;
	       ymid = (y[n_right] + y[n_left]) / 2.;
	       //calculating normals		 
	       nx = yc - ymid;
	       ny = -(xc-xmid);
	       //calculate f+ with Qleft[i] = Q[neqs*n_left+i]
	       n_hat[0] = nx;
	       n_hat[1] = ny;		 
	       f_select[0] = 1; f_select[1] = 0; f_select[2] = 1; f_select[3] = 0;	       
	       calc_van_leer((Q+neqs*n_left), fvl_p, fvl_m, d_fvl_p, d_fvl_m, neqs, gamma, n_hat, f_select);
	       f_select[0] = 0; f_select[1] = 1; f_select[2] = 0; f_select[3] = 1;	       
	       calc_van_leer((Q+neqs*n_right), fvl_p, fvl_m, d_fvl_p, d_fvl_m, neqs, gamma, n_hat, f_select);

	       // accumulating residals to rhs
	       for(j=0; j < neqs; j++)
	       {
		    rhs[neqs*n_left+j] -= (fvl_p[j] + fvl_m[j]);
		    rhs[neqs*n_right+j] += (fvl_p[j] + fvl_m[j]);
	       }
	       //accumulating diagonal
	       accumulate_to_diag(A, d_fvl_p, n_left, neqs, iau);

	       //accumulating off-diagonal
	       accumulate_to_ij(A, d_fvl_m, n_left, n_right, neqs, ia, ja);
	       

	  }
	  
     }

     //clean-up
     free(n_hat);
     free(f_select);     
     free(fvl_p);
     free(fvl_m);
     free(d_fvl_p);
     free(d_fvl_m);
     free(Q_edge);
     free(fw);
     free(dfw);

     //completed successfully
     return 0;
}

//implemets euler implicit scheme in Ax = b = rhs form
int Axb_euler_implicit(double *Q, double *Q_inf, double gamma, double CFL, int ITR_MAX, int itr_per_msg, double *x, double *y, int *bn_nodes, int nn, int neqs, int nt, int **tri_conn, int nnz, int *ia, int *ja,  int *iau, double *A, double *rhs)
{
  int i_start, i_end;
  double *xn = (double *)malloc( neqs * nn * sizeof(double) );
  double *xn1 = (double *)malloc( neqs * nn * sizeof(double) );
  double *x_star = (double *)malloc( neqs * nn * sizeof(double) );
  //locals
  int i,j,k;
  int ITR = 0;
  double *int_uplusc_dl = (double *)calloc(nn , sizeof(double) );

  // main iteration loop
  for( ITR = 1; ITR <= ITR_MAX; ITR++)
    {

      //fill A , rhs
      fill_A_b(Q, Q_inf, gamma, x, y, bn_nodes, nn, neqs, nt, tri_conn, nnz, ia, ja, iau, A, rhs);

      // calculating line integral int( (|u_bar| + c) dl ) 
      calc_int_uplusc_dl( Q, gamma, neqs, nn, x, y, nt, tri_conn, bn_nodes, int_uplusc_dl);

      // add Ai/dt to the diagonal
      for ( i = 0; i < nn; i++)
	{
	  k = iau[i];
	  for( j = 0; j < neqs; j++)
	    A[k*neqs*neqs + j*neqs+j] += int_uplusc_dl[i]/CFL;
	}

      // impose far field BCS DQ = [0]
      for (i = 0; i < nn; i++)
	if( bn_nodes[i] == 1 ) //far field
	  {
	    i_start = ia[i];
	    i_end = ia[i+1]-1;
	    for( j = i_start; j <= i_end; j++)
	      freez((A+ j*neqs*neqs) , neqs, neqs);

	    identity((A+ iau[i]*neqs*neqs), neqs);
	    freez((rhs+i*neqs) , neqs, 1);		    
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
  free(int_uplusc_dl);
  free(xn);
  free(xn1);
  free(x_star);



  //completed successfully!
  return 0;

}
