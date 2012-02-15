#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util2d.h"
//#include "flux.h"
#include "grid_reader.h"
#include "residuals.h"


//the driver routine for fluxes
int main(int argc, char *argv[])
{
     int i,j;
     const int neqs = 4;
     if(argc != 4)
     {
	  printf("\n not enough input arguments! \n syntax:\n $>./unst2d input_mesh_file Mach alpha\n exit ...\n");
	  exit(0);
     }
     double gamma = 1.4;
     double M_inf = atof(argv[2]);
     double alpha = atof(argv[3])*M_PI/180.;
     /* double nx = -1.; */
     /* double ny = 2.0; */
     double *x, *y;
     int nn, nt, **tri_conn, nb, *nbs, ***bs;
     /* printf("\n ---------------- Summary -------------------\n"); */
     /* printf("M_inf = %e, alpha= %e (RAD), nx=%e, ny=%e\n",M_inf, alpha, nx, ny); */

     double *n_hat = (double *)calloc(2 , sizeof(double));
     /* int *f_select = (int *)calloc(4 , sizeof(int)); */

     //allocating vector of conservative variables and Van Leer flux vector 
     double *Q_inf = (double *)calloc( neqs , sizeof(double));
     /* double *fvl_p = (double *)calloc( neqs , sizeof(double)); */
     /* double *fvl_m = (double *)calloc( neqs , sizeof(double)); */

     /* //allocating Van Leer flux vector jacobians d_+-  */
     /* double **d_fvl_p = (double **)calloc( neqs , sizeof(double *)); */
     /* double **d_fvl_m = (double **)calloc( neqs , sizeof(double *)); */

     /* for( i = 0; i < neqs; i++) */
     /* { */
     /* 	  d_fvl_p[i] = (double *)calloc( neqs , sizeof(double)); */
     /* 	  d_fvl_m[i] = (double *)calloc( neqs , sizeof(double)); */
     /* } */

     //initializing Q
     Q_inf[0] = 1.;
     Q_inf[1] = M_inf * cos(alpha);
     Q_inf[2] = M_inf * sin(alpha);
     Q_inf[3] = 1./(gamma * (gamma-1.)) + .5 * M_inf*M_inf;
     //initializing n_hat
     /* n_hat[0] = nx; */
     /* n_hat[1] = ny; */

     //reading the input mesh
     read_mesh_file(argv[1], &x, &y, &nn, &nt, &tri_conn, &nb, &nbs, &bs);

     double *Q = (double *)malloc( neqs * nn * sizeof(double));
     double *R = (double *)malloc( neqs * nn * sizeof(double));

     for (i = 0; i < nn; i++)
	  for (j = 0; j < neqs; j++)
	       Q[i*neqs + j] = Q_inf[j];
     
     /* //calculating fluxes */
     /* f_select[0] = 1;      */
     /* f_select[1] = 1;      */
     /* f_select[2] = 1;      */
     /* f_select[3] = 1; */
     /* calc_van_leer(Q, fvl_p, fvl_m, d_fvl_p, d_fvl_m, neqs, gamma, n_hat, f_select); */
     int *bn_nodes = NULL;     
     tag_bn_nodes(nn, nb, nbs, bs, &bn_nodes);
     
     calc_residuals( Q, Q_inf, gamma, nn, neqs, x, y, nt, tri_conn, bn_nodes, R);

     //showing the matrices
     print_array("residuals",R, neqs*nn);
     printf("\n\n the max(abs(res[j])) = %e\n\n", max_abs_array(R, (neqs*nn)));
     /* print_array("fvl_m",fvl_m, neqs); */

     /* print_matrix("dfplus", d_fvl_p, neqs, neqs); */
     /* print_matrix("dfmin", d_fvl_m, neqs, neqs); */


     //completed successfully!
     return 0;
}
