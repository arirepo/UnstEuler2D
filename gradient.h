#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include <vector>

using namespace std;

int grad_2nd_order(double *Q, vector<int> *p_to_e, int neqs, double *area, double *x, double *y, int *bn_nodes, int nn, int **tri_conn, double *grad);

void test_grad_2nd_order(double *Q, int neqs, double *area, double *x, double *y, int *bn_nodes, int nn, int ntri, int **tri_conn);

#endif
