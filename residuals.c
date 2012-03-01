#include <stdio.h>
#include <stdlib.h>
#include "residuals.h"
#include "flux.h"

//tags the boundary nodes with 1 and interior node with 0.
//so for bn_nodes[i=0..(nn-1)] = 0 --> interior. bn_nodes[i=0..(nn-1)] = 1 --> boundary.     
int tag_bn_nodes(int nn, int nb, int *nbs, int ***bs, int **bn_nodes)
{
     //local vars
     int b, i;
     (*bn_nodes) = (int *)calloc(nn , sizeof(int) );
     
     for (b=0; b < nb; b++)
	  for (i=0; i < nbs[b]; i++)
	  {
	       (*bn_nodes)[bs[b][i][0]] = 1; //tag the first point of boundary
	       (*bn_nodes)[bs[b][i][1]] = 1; //tag the second point of boundary
		     
	  }
     //completed successfully
     return 0;
}

     
//Computes residuals based on median duals
// -------------------- INPUT -----------------------
// Q[neqs*n+i] is the vector of conservative variables [Q0, Qi, ..., Qneqs-1] for node 'n' starting from 0.
// Q_inf[0..3] is the vector of conservative variables at far field. 
// nn is the number of nodes in the grid
// neqs is the number of equations in the system
// x[] and y[] are the coordinates of each node
// nt is the number of triangles in the grid
// tri_conn[tri = i][0...2] shows the nodes that are connected in triangle tri = i according to the predefined counterclockwise winding.
// bn_nodes[node = i] is 1 if the node i (zero based) is a boundary node and zero if it is the interior node. This map should be computed before.
// -------------------- OUTPUT -----------------------
// R[neqs*n+i] is the vector of residuals [R0, Ri, ..., Rneqs-1] for node 'n' starting from 0.

int calc_residuals( double *Q, double *Q_inf, double gamma, int nn, int neqs, double *x, double *y, int nt, int **tri_conn, int *bn_nodes, double *R)
{

     //local vars
     int i, j, t;
     int n_right, n_left;
     double xc, yc, xmid, ymid;
     double nx, ny;

     double *n_hat = (double *)calloc(2 , sizeof(double));
     int *f_select = (int *)calloc(4 , sizeof(int));
     
     double *fvl_p = (double *)calloc( neqs , sizeof(double));
     double *fvl_m = (double *)calloc( neqs , sizeof(double));
     double *Q_edge = (double *)calloc( neqs , sizeof(double));
  
     //allocating Van Leer flux vector jacobians d_+- 
     double **d_fvl_p = (double **)calloc( neqs , sizeof(double *));
     double **d_fvl_m = (double **)calloc( neqs , sizeof(double *));

     for( i = 0; i < neqs; i++)
     {
	  d_fvl_p[i] = (double *)calloc( neqs , sizeof(double));
	  d_fvl_m[i] = (double *)calloc( neqs , sizeof(double));
     }
     //resseting the residuals vector
     for( i = 0; i < nn; i++)	  
	  for(j=0; j < neqs; j++)
	       R[neqs*i+j] = 0.;
     
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

	       if(bn_nodes[n_left] && bn_nodes[n_right]) //both nodes are on the boundaries then this is simply a boundary edge!
	       {
		    nx = 0.5*(y[n_right] - y[n_left]);
		    ny = -0.5*(x[n_right]-x[n_left]);
		    n_hat[0] = nx;
		    n_hat[1] = ny;

		    //for left boundary node
		    for( j = 0; j < neqs; j++)
			 Q_edge[j] = 0.75*Q[neqs*n_left + j] + 0.25*Q[neqs*n_right + j];  		    
		    f_select[0] = 1; f_select[1] = 0; f_select[2] = 0; f_select[3] = 0;	       
		    calc_van_leer(Q_edge, fvl_p, fvl_m, d_fvl_p, d_fvl_m, neqs, gamma, n_hat, f_select);
		    f_select[0] = 0; f_select[1] = 1; f_select[2] = 0; f_select[3] = 0;	       
		    calc_van_leer(Q_inf, fvl_p, fvl_m, d_fvl_p, d_fvl_m, neqs, gamma, n_hat, f_select);
		    for(j=0; j < neqs; j++)
		    {
			 R[neqs*n_left+j] += (fvl_p[j] + fvl_m[j]);
		    }
		    
		    //for right boundary node
		    for( j = 0; j < neqs; j++)
			 Q_edge[j] = 0.75*Q[neqs*n_right + j] + 0.25*Q[neqs*n_left + j];  		    
		    f_select[0] = 1; f_select[1] = 0; f_select[2] = 0; f_select[3] = 0;	       
		    calc_van_leer(Q_edge, fvl_p, fvl_m, d_fvl_p, d_fvl_m, neqs, gamma, n_hat, f_select);
		    f_select[0] = 0; f_select[1] = 1; f_select[2] = 0; f_select[3] = 0;	       
		    calc_van_leer(Q_inf, fvl_p, fvl_m, d_fvl_p, d_fvl_m, neqs, gamma, n_hat, f_select);
		    for(j=0; j < neqs; j++)
		    {
			 R[neqs*n_right+j] += (fvl_p[j] + fvl_m[j]);
		    }
		    
		    
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
	       f_select[0] = 1; f_select[1] = 0; f_select[2] = 0; f_select[3] = 0;	       
	       calc_van_leer((Q+neqs*n_left), fvl_p, fvl_m, d_fvl_p, d_fvl_m, neqs, gamma, n_hat, f_select);
	       f_select[0] = 0; f_select[1] = 1; f_select[2] = 0; f_select[3] = 0;	       
	       calc_van_leer((Q+neqs*n_right), fvl_p, fvl_m, d_fvl_p, d_fvl_m, neqs, gamma, n_hat, f_select);
	       for(j=0; j < neqs; j++)
	       {
		    R[neqs*n_left+j] += (fvl_p[j] + fvl_m[j]);
		    R[neqs*n_right+j] += (-(fvl_p[j] + fvl_m[j]));
	       }
	       

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
     
     //completed successfully
     return 0;
}

//tests if boundary normals all cancel out or not
int test_bn_of_grid( int nn, double *x, double *y, int nt, int **tri_conn, int *bn_nodes, double *sum_nx, double *sum_ny)
{

  //local vars
  int i, j, t;
  int n_right, n_left;
  double xc, yc, xmid, ymid;
  double nx, ny;

  //resetting 
  *sum_nx = 0.;
  *sum_ny = 0.;

  //main loop : loop over all triangles
  for (t = 0; t < nt; t++)
    for(i = 0; i < 3; i++) //loop over vertices of each triangle
      {
	//determine left and right
	n_left = i;	       
	n_right = (i<2)?(i+1):0;

        n_left = tri_conn[t][n_left];
        n_right = tri_conn[t][n_right];

	if(bn_nodes[n_left] && bn_nodes[n_right]) //both nodes are on the boundaries then this is simply a boundary edge!
	  {
	    nx = y[n_right] - y[n_left];
	    ny = -(x[n_right]-x[n_left]);
	    printf("added %lf", nx);
	    printf("added %lf", ny);

	    (*sum_nx) += nx;
	    (*sum_ny) += ny;
	  }
      }

 

  //completed successfully!
  return 0;
}
