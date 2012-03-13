#ifndef _GAUSS_SEIDEL_VALID_H_
#define _GAUSS_SEIDEL_VALID_H_

#define GS_RES_EPS 1.0e-16

int read_A_b_from_file(const char *infile, int *nnodes, int *nnz, int **ia, int **ja,  int **iau, double **A, double **rhs, int neqs);

int gauss_seidel_solve_pivoting(int nnodes, int nnz, int *ia, int *ja,  int *iau, double *A, double *rhs, int neqs, double *x_star, double *xn1, double *xn, short init);

inline void neg_matrix_vec__mult(double *A, double *x,double *x_star,int neqs);

#endif
