#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util2d.h"
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
 
     // check if normal cansel out
     double sum_nx=0., sum_ny=0.;
     test_bn_of_grid( nn, x, y, nt, tri_conn, bn_nodes, &sum_nx, &sum_ny);
    
     printf("sum of nx normal is sum_nx = %19.19e and sum_ny = %19.19e" , sum_nx, sum_ny);
     calc_residuals( Q, Q_inf, gamma, nn, neqs, x, y, nt, tri_conn, bn_nodes, R);

     //showing the matrices
     //print_array("residuals",R, neqs*nn);
     printf("\n\n the max(abs(res[j])) = %e\n\n", max_abs_array(R, (neqs*nn)));
     /* print_array("fvl_m",fvl_m, neqs); */

     /* print_matrix("dfplus", d_fvl_p, neqs, neqs); */
     /* print_matrix("dfmin", d_fvl_m, neqs, neqs); */

//Testing Ariplot
     //1- fill the Q for sample
     for ( i = 0; i < nn; i++)
	  for (j = 0; j < neqs; j++)
	       Q[i*neqs + j] = R[i*neqs + j];
     
     PLT_SPEC samp_plt;
     sprintf(samp_plt.title, "test_contours!");
     sprintf(samp_plt.xlabel, "x");
     sprintf(samp_plt.ylabel, "y");
     samp_plt.xmin = -max_abs_array(x, nn);
     samp_plt.xmax = max_abs_array(x, nn);
     samp_plt.ymin = -max_abs_array(y, nn);
     samp_plt.ymax = max_abs_array(y, nn);
     sprintf(samp_plt.OUTPUT,"display");
     sprintf(samp_plt.pltype, "ColorTri");
     
     write_unst_grd_sol(argv[1], x, y, Q, neqs, nn, nt, tri_conn, &samp_plt);
     
     //completed successfully!
     return 0;
}
