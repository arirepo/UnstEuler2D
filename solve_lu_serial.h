#ifndef _SOLVE_LU_SERIAL_H_
#define _SOLVE_LU_SERIAL_H_

int solve_lu_serial(int *P, double *A, double *x, double *b, int n);
void serial_matrix_mult(int *A, double *B,double *C,int m, int n , int l);

#endif
