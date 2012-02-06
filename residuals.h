#ifndef _RESIDUALS_H_
#define _RESIDUALS_H_

int tag_bn_nodes(int nn, int nb, int *nbs, int ***bs, int **bn_nodes);

int calc_residuals( double *Q, double *Q_inf, double gamma, int nn, int neqs, double *x, double *y, int nt, int **tri_conn, int *bn_nodes, double *R);

#endif
