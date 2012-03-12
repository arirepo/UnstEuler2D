#ifndef _FLUX_H_
#define _FLUX_H_

int calc_van_leer(double *Q, double *fvl_p, double *fvl_m, double **d_fvl_p, double **d_fvl_m, int neqs, double gamma, double *n_hat, int *f_select);

//computes wall flux and Jacobians 
int calc_wall_flux(double *Q, double *fw, double **d_fw, int neqs, double gamma, double *n_hat);

#endif
