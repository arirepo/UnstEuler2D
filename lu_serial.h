#ifndef _LU_SERIAL_H_
#define LU_SERIAL_H_

int lu_serial(double *A, int *P, int n);
inline int find_max(double *A, int *p, int k, int n, int *max_i);
inline int pivot_p(int *p, int i, int j);
void print_matrix_double(const char *name, double *inmat, int n1, int n2);
void print_matrix_int(const char *name, int *inmat, int n1, int n2);
void print_array_double(const char *name, double *inmat, int n1);

#endif

