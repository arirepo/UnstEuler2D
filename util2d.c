
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

//calculates areas per element for 2d unstructure node centered elements (formed by duals)
// 
// output : area[i = 0 ... (nn -1)] = area for element containing node i. 
int cal_areas(int nn, double *x, double *y, int nt, int **tri_conn, double *area)
{
     //local vars
     int i, t;
     int n_right, n_left;
     double xc, yc, xmid, ymid;
     double Ax, Ay, Bx, By;
     double tmp_area = 0.;

     //first resetting area array to zero
     for ( i = 0; i < nn ; i++)
       area[i] = 0.; 
     
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

	       //calculating the center of that edge
	       xmid = (x[n_right] + x[n_left]) / 2.;
	       ymid = (y[n_right] + y[n_left]) / 2.;

	       //now calculate areas for left and right duals and contribute to the left and right nodes.
	       // FIRST : LEFT dual
	       Ax = xmid - x[n_left];
	       Ay = ymid - y[n_left];
	       Bx = xc - x[n_left];
	       By = yc - y[n_left];
	       tmp_area = .5 * fabs(Ax * By - Bx * Ay);
	       area[n_left] += tmp_area;
	       // SECOND : RIGHT dual
	       Ax = xmid - x[n_right];
	       Ay = ymid - y[n_right];
	       Bx = xc - x[n_right];
	       By = yc - y[n_right];
	       tmp_area = .5 * fabs(Ax * By - Bx * Ay);
	       area[n_right] += tmp_area;

	  } //end of loop over vertices
     } //end of loop over all triangles

     //completed successfully
     return 0;

}

//calculates total area of the region occupied by unstructured grid.
// is useful for debugging and validation purposes
// output : total_area = total area of the region
int cal_total_area(int nn, double *x, double *y, int nt, int **tri_conn, double *total_area)
{
  //local vars
  int  t;
  double Ax, Ay, Bx, By;
  double tmp_area = 0.;

  //first resetting total area variable to zero
  *total_area = 0.;
     
  //main loop : loop over all triangles
  for (t = 0; t < nt; t++)
    {

      Ax = x[tri_conn[t][1]] - x[tri_conn[t][0]];
      Ay = y[tri_conn[t][1]] - y[tri_conn[t][0]];
      Bx = x[tri_conn[t][2]] - x[tri_conn[t][0]];
      By = y[tri_conn[t][2]] - y[tri_conn[t][0]];
      tmp_area = .5 * fabs(Ax * By - Bx * Ay);
      (*total_area) += tmp_area;

    } //end of loop over all triangles

  //completed successfully
  return 0;

}

double max_abs_R(double *R, int k, int neqs, int nn)
{
     int i;
     double max = fabs(R[neqs*0+k]);
     for( i = 1; i < nn; i++)
       if(fabs(R[neqs*i+k]) > max)
	       max = fabs(R[neqs*i+k]);

     return max;
     //done!
}

// finds the maximum of an array of integers.
int max_array_int(int *inmat, int n1)
{
     int i;
     int max = inmat[0];
     for( i = 1; i < n1; i++)
	  if(inmat[i] > max)
	       max = inmat[i];

     return max;
     //done!
}

void print_1d_matrix(const char *name, double *inmat, int n1, int n2)
{
     int i,j;
     printf("\n contents of %s : \n", name);
     for( i = 0; i < n1; i++)
     {
	  for( j = 0; j < n2; j++)
	       printf("\t%f,", inmat[i*n2+j]);
	  printf("\n"); 

     }
     //done!
}
