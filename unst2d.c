#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util2d.h"
#include "flux.h"



//the driver routine for fluxes
int main(int argc, char *argv[])
{
     int i;
     const int neqs = 4;
     if(argc != 5)
     {
	  printf("not enough input arguments! \n syntax:\n $>./vanleer Mach alpha nx ny\n exit ...\n");
	  exit(0);
     }
     double gamma = 1.4;
     double M_inf = atof(argv[1]);
     double alpha = atof(argv[2])*M_PI/180.;
     double nx = atof(argv[3]);
     double ny = atof(argv[4]);
     printf("\n ---------------- Summary -------------------\n");
     printf("M_inf = %e, alpha= %e (RAD), nx=%e, ny=%e\n",M_inf, alpha, nx, ny);

     double *n_hat = (double *)calloc(2 , sizeof(double));

     //allocating vector of conservative variables and Van Leer flux vector 
     double *Q = (double *)calloc( neqs , sizeof(double));
     double *fvl_p = (double *)calloc( neqs , sizeof(double));
     double *fvl_m = (double *)calloc( neqs , sizeof(double));

     //allocating Van Leer flux vector jacobians d_+- 
     double **d_fvl_p = (double **)calloc( neqs , sizeof(double *));
     double **d_fvl_m = (double **)calloc( neqs , sizeof(double *));

     for( i = 0; i < neqs; i++)
     {
	  d_fvl_p[i] = (double *)calloc( neqs , sizeof(double));
	  d_fvl_m[i] = (double *)calloc( neqs , sizeof(double));
     }

     //initializing Q
     Q[0] = 1.;
     Q[1] = M_inf * cos(alpha);
     Q[2] = M_inf * sin(alpha);
     Q[3] = 1./(gamma * (gamma-1.)) + .5 * M_inf*M_inf;
     //initializing n_hat
     n_hat[0] = nx;
     n_hat[1] = ny;


     calc_van_leer(Q, fvl_p, fvl_m, d_fvl_p, d_fvl_m, neqs, gamma, n_hat);

     //showing the matrices
     print_array("fvl_p",fvl_p, neqs);
     print_array("fvl_m",fvl_m, neqs);

     print_matrix("dfplus", d_fvl_p, neqs, neqs);
     print_matrix("dfmin", d_fvl_m, neqs, neqs);


     //completed successfully!
     return 0;
}
