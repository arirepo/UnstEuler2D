#ifndef _UTIL2D_H_
#define _UTIL2D_H_

void print_matrix(const char *name, double **inmat, int n1, int n2);
void print_array(const char *name, double *inmat, int n1);
double max_abs_array(double *inmat, int n1);
int cal_areas(int nn, double *x, double *y, int nt, int **tri_conn, double *area);
int cal_total_area(int nn, double *x, double *y, int nt, int **tri_conn, double *total_area);
double max_abs_R(double *R, int k, int neqs, int nn);

int max_array_int(int *inmat, int n1);
void print_1d_matrix(const char *name, double *inmat, int n1, int n2);

#endif
