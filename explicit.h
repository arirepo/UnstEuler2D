#ifndef _EXPLICIT_H_
#define _EXPLICIT_H_

int efficient_euler_explicit(double *Q, double *Q_inf, double gamma, double CFL, int ITR_MAX, int itr_per_msg, int nn, int neqs, double *x, double *y, int nt, int **tri_conn, int *bn_nodes);


#endif
