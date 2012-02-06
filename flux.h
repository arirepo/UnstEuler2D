#ifndef _FLUX_H_
#define _FLUX_H_

int calc_van_leer(double *Q, double *fvl_p, double *fvl_m, double **d_fvl_p, double **d_fvl_m, int neqs, double gamma, double *n_hat);

#endif
