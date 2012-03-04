#ifndef _IMPLICIT_H_
#define _IMPLICIT_H_
#include <vector>
using namespace std;


//allocates sparse matrices A and b required for left and right hand sides of implicit temporal discretization.
int alloc_A_b( int nn, int neqs, int nt, int **tri_conn);





#endif
