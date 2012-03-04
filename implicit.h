#ifndef _IMPLICIT_H_
#define _IMPLICIT_H_
using namespace std;


//allocates sparse matrices A and b required for left and right hand sides of implicit temporal discretization.
int alloc_A_b( int nn, int neqs, int nt, int **tri_conn, int *nnz, int **ia, int **ja,  int **iau, double **A, double **rhs);

//implemets euler explicit scheme in Ax = b = rhs form
int Axb_euler_explicit(double *Q, double *Q_inf, double gamma, double CFL, int ITR_MAX, int itr_per_msg, double *x, double *y, int *bn_nodes, int nn, int neqs, int nt, int **tri_conn, int nnz, int *ia, int *ja,  int *iau, double *A);


#endif
