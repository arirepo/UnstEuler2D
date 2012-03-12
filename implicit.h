#ifndef _IMPLICIT_H_
#define _IMPLICIT_H_
using namespace std;


//allocates sparse matrices A and b required for left and right hand sides of implicit temporal discretization.
int alloc_A_b( int nn, int neqs, int nt, int **tri_conn, int *nnz, int **ia, int **ja,  int **iau, double **A, double **rhs);

//implemets euler explicit scheme in Ax = b = rhs form
int Axb_euler_explicit(double *Q, double *Q_inf, double gamma, double CFL, int ITR_MAX, int itr_per_msg, double *x, double *y, int *bn_nodes, int nn, int neqs, int nt, int **tri_conn, int nnz, int *ia, int *ja,  int *iau, double *A, double *rhs);

//accumulate Df sub-block to the diagonal of matrix A located at (row, row)
inline int accumulate_to_diag(double *A, double **Df, int row, int neqs, int *iau);

//accumulate Df sub-block to the locaation (i,j) of matrix A.
inline int accumulate_to_ij(double *A, double **Df, int i, int j, int neqs, int *ia, int *ja);

// first reset A and b
// fills  martices A and b with Jacobians and residuals respectively
int fill_A_b(double *Q, double *Q_inf, double gamma, double *x, double *y, int *bn_nodes, int nn, int neqs, int nt, int **tri_conn, int nnz, int *ia, int *ja,  int *iau, double *A, double *rhs);


//implemets euler implicit scheme in Ax = b = rhs form
int Axb_euler_implicit(double *Q, double *Q_inf, double gamma, double CFL, int ITR_MAX, int itr_per_msg, double *x, double *y, int *bn_nodes, int nn, int neqs, int nt, int **tri_conn, int nnz, int *ia, int *ja,  int *iau, double *A, double *rhs);

inline void identity(double *S, int neqs);
inline void freez(double *S, int n1, int n2 );


#endif
