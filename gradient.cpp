#include <stdio.h>
#include "gradient.h"
#include <vector>
#include "maps.h"

using namespace std;


//calculates gradients using Divergence theorem.
// stores grad_x and grad_y in the following form:
// grad[i*neqs*2+k*2 + 0] = grad_x for node i equation number k
// grad[i*neqs*2+k*2 + 1] = grad_y for node i equation number k

int grad_2nd_order(double *Q, vector<int> *p_to_e, int neqs, double *area, double *x, double *y, int *bn_nodes, int nn, int **tri_conn, double *grad)
{

  //locals
  int i, j, k, t;
  double nx, ny;
  int right, left;
  int vtri = 3;
  int p1 = 0, p2 = 0;

  //resetting the gradient
  for( i = 0; i < nn; i++)
    for( j = 0; j < neqs; j++)
      for( k = 0; k < 2; k++)
	grad[i*neqs*2+j*2+k] = 0.;

  //main loop over nodes this time
  for ( i = 0; i < nn ; i++)
    {

      for ( j = 0; j < p_to_e[i].size() ; j++) //loop over elements surrounding point
	{
	  t = p_to_e[i][j]; //get that element

	  //finding the node in triangle
	  for ( k = 0; k < 3 ; k++)
	    if( i == tri_conn[t][k] ) //I got it!
	      {
		left = ((k-1)%vtri+vtri)%vtri;
		right = (k+1)%vtri;
		//converting local left and right to global left and right
		p2 = tri_conn[t][left]; //this is left 
		p1 = tri_conn[t][right]; //this is right
	      }

	  // normal to that triangle
	  nx = y[p2] - y[p1];
	  ny = -(x[p2] - x[p1]);

	  //contributing to gradient
	  for ( k = 0; k < neqs ; k++)
	    {
	      grad[i*neqs*2+k*2 + 0] += .5*(Q[p1*neqs + k] + Q[p2*neqs + k])*nx;
	      grad[i*neqs*2+k*2 + 1] += .5*(Q[p1*neqs + k] + Q[p2*neqs + k])*ny;
	    }

	  if( (bn_nodes[i] > 0) && (bn_nodes[p1] > 0) )  //then they are on the boundaries 
	    {
	      //normal to the right boundary edge
	      nx = y[p1] - y[i];
	      ny = -(x[p1] - x[i]);

	      for ( k = 0; k < neqs ; k++)
		{
		  grad[i*neqs*2+k*2 + 0] += .5*(Q[i*neqs + k] + Q[p1*neqs + k])*nx;
		  grad[i*neqs*2+k*2 + 1] += .5*(Q[i*neqs + k] + Q[p1*neqs + k])*ny;
		}
	    }
	  if( (bn_nodes[i] > 0) && (bn_nodes[p2] > 0) )  //then they are on the boundaries 
	    {
	      //normal to the right boundary edge

	      nx = y[i] - y[p2];
	      ny = -(x[i] - x[p2]);

	      for ( k = 0; k < neqs ; k++)
		{
		  grad[i*neqs*2+k*2 + 0] += .5*(Q[i*neqs + k] + Q[p2*neqs + k])*nx;
		  grad[i*neqs*2+k*2 + 1] += .5*(Q[i*neqs + k] + Q[p2*neqs + k])*ny;
		}


	    }
     
      
	} // end of loop over surrounding triangles for each node

      // dividing by areas
      for ( k = 0; k < neqs ; k++)
	{
	  grad[i*neqs*2+k*2 + 0] /= (3.*area[i]);
	  grad[i*neqs*2+k*2 + 1] /= (3.*area[i]);
	}
    
    }// end of loop for all nodes. gradients will be calculated after this 
  //for all nodes.
   


  //completed successfully!
  return 0;

}

void test_grad_2nd_order(double *Q, int neqs, double *area, double *x, double *y, int *bn_nodes, int nn, int ntri, int **tri_conn)
{

  //locals
  int i,j;

  //allocating p_to_e map ...
  vector<int> *p_to_e = (vector<int> *)calloc( nn , sizeof(vector<int>));
  double *grad = (double *)calloc( nn*neqs*2 , sizeof(double));


  //make p2e map 
  create_p_to_e(ntri, tri_conn, p_to_e);

  // initialize the Q
  for(i=0; i < nn; i++)
    for(j=0; j < neqs; j++)
      Q[i*neqs + j] = (double)j*x[i] + ((double)j+1.)*y[i];
  // find gradients
  grad_2nd_order(Q, p_to_e, neqs, area, x, y, bn_nodes, nn, tri_conn, grad);

  //print gradients
  for(i=0; i < nn; i++)
    for(j=0; j < neqs; j++)
      printf("node %d : grad_x%d = %e, grad_y%d = %e\n", i, j, grad[i*neqs*2+j*2+0], j, grad[i*neqs*2+j*2+1]);


}

