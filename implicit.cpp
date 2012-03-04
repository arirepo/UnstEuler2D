#include <stdio.h>
#include "implicit.h"
#include "maps.h"
#include <vector>

using namespace std;

//allocates sparse matrices A and b required for left and right hand sides of implicit temporal discretization.
int alloc_A_b( int nn, int neqs, int nt, int **tri_conn)
{
  //local indices
  int i,j; 

  //allocating p2p map ...
  vector<int> *p2p = (vector<int> *)calloc( nn , sizeof(vector<int>));
  // creating p2p map ...
  create_p2p(nt, tri_conn, nn, p2p);


  //debugging
  for( i = 0; i < nn; i++)
    { 
      printf("\n p2p[%d]: ", i);
      fflush(stdout);
      for ( j = 0; j < p2p[i].size(); j++)
	printf("%d,", p2p[i][j]);
    }
  printf("\n");


  //completed successfully!
  return 0;

}
