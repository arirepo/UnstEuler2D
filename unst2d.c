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
     double *x, *y;
     int nn, nt, **tri_conn, nb, *nbs, ***bs;
     double *Q_inf = (double *)calloc( neqs , sizeof(double));

     //initializing Q
     Q_inf[0] = 1.;
     Q_inf[1] = M_inf * cos(alpha);
     Q_inf[2] = M_inf * sin(alpha);
     Q_inf[3] = 1./(gamma * (gamma-1.)) + .5 * M_inf*M_inf;

     //reading the input mesh
     read_mesh_file(argv[1], &x, &y, &nn, &nt, &tri_conn, &nb, &nbs, &bs);

     double *Q = (double *)malloc( neqs * nn * sizeof(double));
     double *R = (double *)malloc( neqs * nn * sizeof(double));

     for (i = 0; i < nn; i++)
	  for (j = 0; j < neqs; j++)
	    Q[i*neqs + j] = (Q_inf[j]+.01);

     //tag the boundary nodes     
     int *bn_nodes = NULL;     
     tag_bn_nodes(nn, nb, nbs, bs, &bn_nodes);
 
     // check if normal cansel out
     double sum_nx=0., sum_ny=0.;
     test_bn_of_grid( nn, x, y, nt, tri_conn, bn_nodes, &sum_nx, &sum_ny);
    
     printf("\nsum of nx normal is sum_nx = %19.19e and sum_ny = %19.19e\n" , sum_nx, sum_ny);

     //calculating areas
     double *area = (double *)malloc( nn * sizeof(double));
     double sum_indv_areas = 0.;
     double total_area = 0.;
     // calculating area per node-centered element
     cal_areas(nn, x, y, nt, tri_conn, area);
     for ( i = 0 ; i < nn; i++)
       sum_indv_areas += area[i];

     cal_total_area(nn, x, y, nt, tri_conn, &total_area);
     printf("\ntotal area is = %17.17e and sum of area array is %17.17e\n" , total_area, sum_indv_areas);

     int ITR = 0;
     double *int_uplusc_dl = (double *)malloc(nn * sizeof(double) );
     double CFL = 25.2;

     // main iteration loop
     for( ITR = 1; ITR < 1000; ITR++)
       {

	 //finding the residuals
	 calc_residuals( Q, Q_inf, gamma, nn, neqs, x, y, nt, tri_conn, bn_nodes, R);

	 // calculating line integral int( (|u_bar| + c) dl ) 
	 calc_int_uplusc_dl( Q, gamma, neqs, nn, x, y, nt, tri_conn, bn_nodes, int_uplusc_dl);

	 printf("ITR = %d, norm(R) = %17.17e\n", ITR, max_abs_array(R, (neqs*nn)));
	 //updating Q
	 for( i = 0; i < nn; i++)
	   for( j = 0; j < neqs; j++)
	     {
	       //printf("DQ = %e\n", CFL*area[i]/int_uplusc_dl[i]);
	       Q[i*neqs + j] = Q[i*neqs + j] - CFL*area[i]/int_uplusc_dl[i] * R[i*neqs + j];
	     }


       }

     printf("\n\n the max(abs(res[j])) = %e\n\n", max_abs_array(R, (neqs*nn)));
     print_array("Q_inf", Q_inf, 4);

     //Testing Ariplot
     
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
