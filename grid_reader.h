#ifndef _GRID_READER_H_
#define _GRID_READER_H_

int read_mesh_file(char *fname, double **x, double **y, int *nn, int *nt, int ***tri_conn, int *nb, int **nbs, int ****bs);

int write_mesh_gnu_files(char *fname, double *x, double *y, int nn, int nt, int tri[][3], int nb, int *nbs, int ***bs);


#endif
