#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util2d.h"
#include "grid_reader.h"
#include "residuals.h"
#include "explicit.h"
#include "implicit.h"

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
	    Q[i*neqs + j] = Q_inf[j];

     //tag the boundary nodes     
     int *bn_nodes = NULL;     
     tag_bn_nodes(nn, nb, nbs, bs, &bn_nodes);
     /* // plot boundary node numbers */
     /* double *bn_nodes_plot = (double *)malloc(nn *neqs * sizeof(double)); */
     /* for( i = 0; i < nn; i++) */
     /*   for (j = 0 ; j < neqs; j++) */
     /*   bn_nodes_plot[i*neqs+j] = (double)bn_nodes[i]; */
     /* { */
     /*   PLT_SPEC samp_plt; */
     /*   sprintf(samp_plt.title, "boundary_number_plus_one"); */
     /*   sprintf(samp_plt.xlabel, "x"); */
     /*   sprintf(samp_plt.ylabel, "y"); */
     /*   samp_plt.xmin = -max_abs_array(x, nn); */
     /*   samp_plt.xmax = max_abs_array(x, nn); */
     /*   samp_plt.ymin = -max_abs_array(y, nn); */
     /*   samp_plt.ymax = max_abs_array(y, nn); */
     /*   sprintf(samp_plt.OUTPUT,"display"); */
     /*   sprintf(samp_plt.pltype, "Contour"); */
     
     /*   write_unst_grd_sol(argv[1], x, y, bn_nodes_plot, neqs, nn, nt, tri_conn, &samp_plt); */
     /*   //print_array("bn_nodes_plot", bn_nodes_plot, nn*neqs); */
     /*   printf("boundary number visualization finished! exit ...\n"); */
     /*   exit(0); */
     /* } */

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

     //starting to march with matrix independent implementation of euler explicit scheme
     double CFL = .9;
     int ITR_MAX = 15000;
     int itr_per_msg = 1;
     //efficient_euler_explicit(Q, Q_inf, gamma, CFL, ITR_MAX, itr_per_msg, nn, neqs, x, y, nt, tri_conn, bn_nodes);

     // allocating sparse matrices [A] and [b] based on the input grid
     int nnz = 0;
     int *ia = NULL, *ja = NULL, *iau = NULL;
     double *A = NULL, *rhs = NULL;
     alloc_A_b( nn, neqs, nt, tri_conn, &nnz, &ia, &ja,  &iau, &A, &rhs);
     //visualize node numbers in gnuplot 
     xy_tri_gnu_plot("sample_node_number.dat", x, y, tri_conn, nt);

     Axb_euler_explicit(Q, Q_inf, gamma, CFL, ITR_MAX, itr_per_msg, x, y, bn_nodes, nn, neqs, nt, tri_conn, nnz, ia, ja, iau, A);

     //Testing Ariplot     
     PLT_SPEC samp_plt;
     sprintf(samp_plt.title, "contours_e_Minf%1.1f_alpha%1.1f", M_inf, alpha*180./M_PI);
     sprintf(samp_plt.xlabel, "x");
     sprintf(samp_plt.ylabel, "y");
     samp_plt.xmin = -max_abs_array(x, nn);
     samp_plt.xmax = max_abs_array(x, nn);
     samp_plt.ymin = -max_abs_array(y, nn);
     samp_plt.ymax = max_abs_array(y, nn);
     sprintf(samp_plt.OUTPUT,"display");
     sprintf(samp_plt.pltype, "Contour");
     
     write_unst_grd_sol(argv[1], x, y, Q, neqs, nn, nt, tri_conn, &samp_plt);

     //clean - ups 
     free(x);
     free(y);
     free(tri_conn);
     free(nbs);
     free(bs);

     free(Q_inf);
     free(Q);
     free(R);
     free(bn_nodes);
     free(area);
     
     //completed successfully!
     return 0;
}
