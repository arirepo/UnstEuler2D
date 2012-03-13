#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gauss_seidel_valid.h"
#include "lu_serial.h"
#include "solve_lu_serial.h"

int read_A_b_from_file(const char *infile, int *nnodes, int *nnz, int **ia, int **ja,  int **iau, double **A, double **rhs, int neqs)
{
     //local vars
     int i = 0, j = 0, l = 0, m = 0;
     int jstart, jend;
     FILE *fid = NULL;
     const int str_dim = 500;
     char str[str_dim];
     
     //openning the file
     if( (fid = fopen(infile, "r")) == NULL)
     {
      printf ("\ncouldn't read from input file. Error occured at file = %s, line %d. Proceeding to exit ...\n",__FILE__, __LINE__);
      fflush(stdout);
      exit(0);
     }
     
     //starting to read
     fgets(str,str_dim,fid); //skip that line  
     fscanf(fid,"%d\n",nnodes); //number of nodes
     printf("\n nnodes = %d", *nnodes);
     fgets(str,str_dim,fid); //skip that line 
     fscanf(fid,"%d\n",nnz); //number of entries
     printf("\n nnz = %d", *nnz);     
     //allocating required arrays
     //THESE WILL BE RETURNED BACK
     (*ia) = (int *)malloc( (*nnodes + 1) * sizeof(int) );
     (*iau) = (int *)malloc( (*nnodes) * sizeof(int) );
     (*ja) = (int *)malloc( (*nnz) * sizeof(int) );
     (*A) = (double *)malloc( neqs * neqs * (*nnz) * sizeof(double) );
     (*rhs) = (double *)malloc( neqs * (*nnodes) * sizeof(double) );

     fgets(str,str_dim,fid); //skip that line     
     //reading ia
     for(i = 0; i < (*nnodes + 1); i++)
     {
	  fscanf(fid,"%d\n",(*ia + i));
	  (*ia)[i]--; //resetting zero-based!	  
     }
     printf("\n ia[end] = %d", (*ia)[*nnodes]);     
     fgets(str,str_dim,fid); //skip that line
     //reading iau
     for(i = 0; i < (*nnodes) ; i++)
     {
	  fscanf(fid,"%d\n",(*iau + i));
	  (*iau)[i]--; //resetting zero-based!	  
     }
     printf("\n iau[end] = %d", (*iau)[*nnodes-1]);     
     fgets(str,str_dim,fid); //skip that line
     //reading ja
     for(i = 0; i < (*nnz) ; i++)
     {
	  fscanf(fid,"%d\n",(*ja + i));
	  (*ja)[i]--; //resetting zero-based!
     }
     printf("\n ja[end] = %d", (*ja)[*nnz-1]);          
     //reading sparse matrix [A]
     fgets(str,str_dim,fid); //skip that line     
     for( i = 0; i < (*nnodes); i++)
     {
	  jstart = (*ia)[i];
	  jend =  (*ia)[i+1] - 1;
	  for (j = jstart; j <= jend; j++)	 
	       for(l = 0; l < neqs; l++)
	       {
		    for(m = 0; m < (neqs-1); m++)
			 fscanf(fid,"%lf ", (*A + j* neqs * neqs+ l * neqs + m));
		    fscanf(fid,"%lf\n", (*A + j* neqs * neqs+ l * neqs + m));
	       }
     }
     printf("\n A[end,end,end] = %1.16f", *(*A + (j-1)* neqs * neqs+ (l-1) * neqs + m)); 
     fgets(str,str_dim,fid); //skip that line
     for( i = 0; i < (*nnodes); i++)
     {

	  for( m = 0; m < (neqs-1); m++)
	       fscanf(fid,"%lf ",(*rhs + i* neqs + m));

	  fscanf(fid,"%lf\n", (*rhs + i* neqs + m));
     }
     printf("\n rhs[end,end] = %1.16f", *(*rhs + (i-1)* neqs + m)); 

     //closing the input file
     fclose(fid);

     printf("\nMatrices A and rhs were ported successfully!\n");
     
     //completed successfully!
     return 0;
}


//solves Ax = b using Gauss-Seidel algorithm
     // --------------- p s e u d o c o d e -----------------------
     //replace the diagonal matrices with their LU and save in-place
     //initial guess [xn] = [x]_0
     //do {
     //    calculate [x_star] = -[O] [xn]
     //    [x_star] = [rhs]  + [x_star]
     //    [xn+1] = [D]^-1 [x_star] (lu-solve element by element)
     //} while (norm([xn+1] - [xn]) > epsilon)
     //export solution
     //exit
     // ----------------------------------------------------------

//ver 1.00, copyright 2012 A. Ghasemi, ghasemi.arash@gmail.com 
//license : BSD - see file attached -

int gauss_seidel_solve_pivoting(int nnodes, int nnz, int *ia, int *ja,  int *iau, double *A, double *rhs, int neqs, double *x_star, double *xn1, double *xn, short init)
{
     //locals
     int i,j;
     int jstart, jend;
     double GS_res = 0.;
     
     // allocating permutation matrix
     int *P = (int *)malloc( nnodes * neqs * neqs * sizeof(int) );
 
     //replace the diagonal matrices with their LU and save in-place
     //also save permutation matrix band per main diagonal in [P]
     for ( i = 0; i < nnodes; i++)
	  lu_serial((A + iau[i]*neqs*neqs), (P + i*neqs*neqs), neqs);

     //initializing solution [xn] with [rhs]/A(diag)
     if (init)
       for ( i = 0; i < nnodes; i++)
	 for( j = 0; j < neqs; j++)
	   {
	     jstart = iau[i];
	     xn[i*neqs + j] = rhs[i*neqs + j] / A[jstart*neqs*neqs + j*neqs + j];
	   }
     
     do{
	  //resetting [x_star]
	  for ( i = 0; i < nnodes; i++)
	       for( j = 0; j < neqs; j++)
		    x_star[i*neqs + j] = 0.;
	  
	  //calculating [x_star] = -[O] [xn].
	  for( i = 0; i < nnodes; i++)
	  {
	       jstart = ia[i];
	       jend = ia[i+1]-1;
	       for ( j = jstart; j <= jend; j++)
		    if( j == iau[i] ) // this is diagonal block, skip it!
			 continue;
		    else //contribute the off-diagonal elements
			 neg_matrix_vec__mult( (A + j* neqs *neqs) , (xn + ja[j]*neqs) , (x_star + i *neqs) , neqs);
	       
	       for ( j = 0; j < neqs; j++) //[x_star] = [rhs]  + [x_star]
		    x_star[i *neqs + j] += rhs[i * neqs + j];
		    
	  } //matrix multiplication and addition will be finished after this!

	  //solving diagonal, i.e. [xn+1] = [D]^-1 [x_star]
	  for( i = 0; i < nnodes; i++)
	  {
	       j = iau[i];
	       solve_lu_serial((P + i*neqs*neqs), (A + j*neqs*neqs), (xn1 + i*neqs), (x_star + i*neqs), neqs);
	       
	  } //xn+1 will be obtained after this loop
	  //calculating residuals
	  GS_res = 0.; //initial
	  for ( i = 0; i < nnodes; i++)
	       for( j = 0; j < neqs; j++)
		    GS_res += pow((xn1[i*neqs + j] - xn[i*neqs + j]) , 2.);
	  GS_res /= (neqs * nnodes);
	  GS_res = sqrt(GS_res);

	  //updating [xn]
	  for ( i = 0; i < nnodes; i++)
	       for( j = 0; j < neqs; j++)
		    xn[i*neqs + j]  = xn1[i*neqs + j];
	  //printf("\n%e\n", GS_res);	  
     }while( GS_res >= GS_RES_EPS );

     //clean - up
     free(P);
     
     //completed successfully!
     return 0;

}

/*      Performs the NEGATIVE  matrix-vector multipication
	and ACCUMULATES the result to x_star 
                        x_star += (-A*x)
	for A[neqs x neqs] and x[neqs x 1]  		
 */
inline void neg_matrix_vec__mult(double *A, double *x,double *x_star,int neqs)
{
     /* local variables */
     int i,j;
     for ( i= 0; i < neqs ; i++)
	  for ( j= 0; j < neqs ; j++)
	       x_star[i] -= A[i*neqs+j] * x[j];

     /* simply finished .*/
     
}
