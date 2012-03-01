#ifndef _GRID_READER_H_
#define _GRID_READER_H_

typedef struct {
     char font[50];
     double font_size;
     char title[500];
     char xlabel[500];
     char ylabel[500];
     double xmin;
     double xmax;
     double ymin;
     double ymax;
     char OUTPUT[500];
     char pltype[100];
          
} PLT_SPEC;

int read_mesh_file(char *fname, double **x, double **y, int *nn, int *nt, int ***tri_conn, int *nb, int **nbs, int ****bs);

int write_unst_grd_sol(char *fname, double *x, double *y, double *Q, int neqs, int nn, int nt, int **tri, PLT_SPEC *gnplt);



#endif
