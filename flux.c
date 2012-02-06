#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//prototyping
void print_matrix(const char *name, double **inmat, int n1, int n2);
int calc_van_leer(double *Q, double *fvl_p, double *fvl_m, double **d_fvl_p, double **d_fvl_m, int neqs, double gamma, double *n_hat);
void print_array(const char *name, double *inmat, int n1);

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

void print_matrix(const char *name, double **inmat, int n1, int n2)
{
     int i,j;
     printf("\n contents of %s : \n", name);
     for( i = 0; i < n1; i++)
     {
	  for( j = 0; j < n2; j++)
	       printf("%1.16f, ", inmat[i][j]);
	  printf("\n"); 

     }
     //done!
}

void print_array(const char *name, double *inmat, int n1)
{
     int i;
     printf("\n contents of %s : \n", name);
     for( i = 0; i < n1; i++)
	  printf("%1.16f \n", inmat[i]);
     //done!
}

int calc_van_leer(double *Q, double *fvl_p, double *fvl_m, double **d_fvl_p, double **d_fvl_m, int neqs, double gamma, double *n_hat)
{
     int i, j;
     double sign_pm = 1., p2 = 0.;
     double length = 0.;
     double nx, ny;
     double rho, u, v, u_bar, e, P, c;
     double *d_rho_dQi = (double *)calloc(neqs , sizeof(double));
     double *d_u_dQi = (double *)calloc(neqs , sizeof(double));
     double *d_v_dQi = (double *)calloc(neqs , sizeof(double));
     double *d_e_dQi = (double *)calloc(neqs , sizeof(double));
     double *d_P_dQi = (double *)calloc(neqs , sizeof(double));
     double *d_c_dQi = (double *)calloc(neqs , sizeof(double));
     double *d_ubar_dQi = (double *)calloc(neqs , sizeof(double));  
 
     //step 1: calculating primitive variables
     length = sqrt(n_hat[0] * n_hat[0]+ n_hat[1] * n_hat[1]);
     nx = n_hat[0]/length;
     ny = n_hat[1]/length;
     rho = Q[0];
     u = Q[1] / Q[0];
     v = Q[2] / Q[0];
     e = Q[3];
     u_bar = u * nx + v * ny;
     P = (gamma - 1.) * e - .5 * (gamma - 1.) *rho * ( u*u + v*v);
     c = sqrt(gamma * P/ rho);

     //step 2: filling the flux van leer (fvl) vector ASAP
     if(u_bar >= c)
     {
	  fvl_p[0] = length * rho* u_bar;
	  fvl_p[1] = length * (rho * u * u_bar + nx * P);
	  fvl_p[2] = length * (rho * v * u_bar + ny * P);
	  fvl_p[3] = length * (e+P)* u_bar;

	  fvl_m[0] = 0.;
	  fvl_m[1] = 0.;
	  fvl_m[2] = 0.;
	  fvl_m[3] = 0.;

     }
     else if(u_bar <= (-c))
     {
	  fvl_m[0] = length * rho* u_bar;
	  fvl_m[1] = length * (rho * u * u_bar + nx * P);
	  fvl_m[2] = length * (rho * v * u_bar + ny * P);
	  fvl_m[3] = length * (e+P)* u_bar;

	  fvl_p[0] = 0.;
	  fvl_p[1] = 0.;
	  fvl_p[2] = 0.;
	  fvl_p[3] = 0.;

     }
     else
     {
	  fvl_p[0] = length * .250 * rho*c* pow((u_bar/c + 1.) , 2.);
	  fvl_p[1] = fvl_p[0] * (nx/gamma*(-u_bar + 2. * c) + u);
	  fvl_p[2] = fvl_p[0] * (ny/gamma*(-u_bar + 2. * c) + v);
	  fvl_p[3] = fvl_p[0] * ( (-(gamma-1)*u_bar*u_bar + 2.*(gamma-1.)*u_bar*c+2.*c*c)/ (gamma*gamma - 1.) + (u*u + v*v)*.5 );

	  fvl_m[0] = length * -.250 * rho*c* pow((u_bar/c - 1.) , 2.);
	  fvl_m[1] = fvl_m[0] * (nx/gamma*(-u_bar - 2. * c) + u);
	  fvl_m[2] = fvl_m[0] * (ny/gamma*(-u_bar - 2. * c) + v);
	  fvl_m[3] = fvl_m[0] * ( (-(gamma-1)*u_bar*u_bar - 2.*(gamma-1.)*u_bar*c+2.*c*c)/ (gamma*gamma - 1.) + (u*u + v*v)*.5 );

     }
  
     //step 3: calculating derivatives of primitive variables
     d_rho_dQi[0] = 1.;   d_rho_dQi[1] = 0.;   d_rho_dQi[2] = 0.;   d_rho_dQi[3] = 0.;
     d_u_dQi[0] = -Q[1]/(Q[0]*Q[0]); d_u_dQi[1] = 1./Q[0]; d_u_dQi[2] = 0.; d_u_dQi[3] = 0.;
     d_v_dQi[0] = -Q[2]/(Q[0]*Q[0]); d_v_dQi[1] = 0.; d_v_dQi[2] = 1./Q[0]; d_v_dQi[3] = 0.;
     d_e_dQi[0] = 0.; d_e_dQi[1] = 0.; d_e_dQi[2] = 0.; d_e_dQi[3] = 1.;

     //step 4: calculating derivatives of physical variables P, c, u_bar
  
     for( i = 0; i < neqs; i++)
     {
	  d_P_dQi[i] = (gamma-1.) * (d_e_dQi[i] - .5 * d_rho_dQi[i] *(u*u + v*v) - rho * (u * d_u_dQi[i] + v * d_v_dQi[i]));
	  d_c_dQi[i] = gamma/2. * (d_P_dQi[i] * rho - d_rho_dQi[i] * P) / (rho * rho * c);
	  d_ubar_dQi[i] = d_u_dQi[i] * nx + d_v_dQi[i] * ny;  
     }

     //step 5: using chain-rule to finish the job.
     // starting filling row by row ..
     if(u_bar >= c)
     {
	  // d_fvl_p
	  for( j = 0; j < neqs; j++)
	  {
	       d_fvl_p[0][j] = length*(d_rho_dQi[j]*u_bar + rho*d_ubar_dQi[j]);
	       d_fvl_p[1][j] = length * (d_rho_dQi[j] * u * u_bar+ rho * d_u_dQi[j] * u_bar+ rho * u * d_ubar_dQi[j] + nx * d_P_dQi[j]); 
	       d_fvl_p[2][j] = length * (d_rho_dQi[j] * v * u_bar+ rho * d_v_dQi[j] * u_bar+ rho * v * d_ubar_dQi[j] + ny * d_P_dQi[j]); 
	       d_fvl_p[3][j] = length * ( (d_e_dQi[j]+d_P_dQi[j])* u_bar + (e+P)* d_ubar_dQi[j] ); 
	  }

	  // d_fvl_m
	  for( j = 0; j < neqs; j++)
	  {
	       d_fvl_m[0][j] = 0.;
	       d_fvl_m[1][j] = 0.;
	       d_fvl_m[2][j] = 0.;
	       d_fvl_m[3][j] = 0.;
	  }
     }
     else if(u_bar <= (-c))
     {
	  // d_fvl_p
	  for( j = 0; j < neqs; j++)
	  {
	       d_fvl_p[0][j] = 0.;
	       d_fvl_p[1][j] = 0.;
	       d_fvl_p[2][j] = 0.;
	       d_fvl_p[3][j] = 0.;
	  }
	  // d_fvl_m
	  for( j = 0; j < neqs; j++)
	  {
	       d_fvl_m[0][j] = length*(d_rho_dQi[j]*u_bar + rho*d_ubar_dQi[j]);
	       d_fvl_m[1][j] = length * (d_rho_dQi[j] * u * u_bar+ rho * d_u_dQi[j] * u_bar+ rho * u * d_ubar_dQi[j] + nx * d_P_dQi[j]); 
	       d_fvl_m[2][j] = length * (d_rho_dQi[j] * v * u_bar+ rho * d_v_dQi[j] * u_bar+ rho * v * d_ubar_dQi[j] + ny * d_P_dQi[j]); 
	       d_fvl_m[3][j] = length * ( (d_e_dQi[j]+d_P_dQi[j])* u_bar + (e+P)* d_ubar_dQi[j] ); 
	    
	  }       
     }
     else
     {
	  // d_fvl_p
	  sign_pm = 1.0;
	  p2 = pow((u_bar / c + sign_pm * 1.),2.0); 
	  for( j = 0; j < neqs; j++)
	  {
	       d_fvl_p[0][j] = length*sign_pm*.25*(d_rho_dQi[j]*c*p2 + rho * d_c_dQi[j] * p2 + 2.*(rho/c)* (u_bar/c + sign_pm * 1.) * (d_ubar_dQi[j]*c-d_c_dQi[j]*u_bar));
	       d_fvl_p[1][j] = d_fvl_p[0][j] * (nx/gamma * (-u_bar + sign_pm*2.*c)+u) + fvl_p[0] * (nx/gamma*(-d_ubar_dQi[j]+sign_pm* 2.* d_c_dQi[j])+ d_u_dQi[j]); 
	       d_fvl_p[2][j] = d_fvl_p[0][j] * (ny/gamma * (-u_bar + sign_pm*2.*c)+v) + fvl_p[0] * (ny/gamma*(-d_ubar_dQi[j]+sign_pm* 2.* d_c_dQi[j])+ d_v_dQi[j]); 
	       d_fvl_p[3][j] = d_fvl_p[0][j] * ((-(gamma-1.)*u_bar*u_bar+ sign_pm*2.*(gamma-1.)*u_bar*c+2.*c*c)/(gamma*gamma-1.)+.5*(u*u+v*v)) + fvl_p[0] * ((-2.*(gamma-1.)*u_bar*d_ubar_dQi[j]+sign_pm*2.*(gamma-1.)*(d_ubar_dQi[j]*c+u_bar*d_c_dQi[j])+ 4.*c*d_c_dQi[j])/(gamma*gamma-1.)+ u*d_u_dQi[j] + v*d_v_dQi[j]); 
	  }

	  // d_fvl_m
	  sign_pm = -1.0;
	  p2 = pow((u_bar / c + sign_pm * 1.),2.0); 
	  for( j = 0; j < neqs; j++)
	  {
	       d_fvl_m[0][j] = length*sign_pm*.25*(d_rho_dQi[j]*c*p2 + rho * d_c_dQi[j] * p2 + 2.*(rho/c)* (u_bar/c + sign_pm * 1.) * (d_ubar_dQi[j]*c-d_c_dQi[j]*u_bar));
	       d_fvl_m[1][j] = d_fvl_m[0][j] * (nx/gamma * (-u_bar + sign_pm*2.*c)+u) + fvl_m[0] * (nx/gamma*(-d_ubar_dQi[j]+sign_pm* 2.* d_c_dQi[j])+ d_u_dQi[j]); 
	       d_fvl_m[2][j] = d_fvl_m[0][j] * (ny/gamma * (-u_bar + sign_pm*2.*c)+v) + fvl_m[0] * (ny/gamma*(-d_ubar_dQi[j]+sign_pm* 2.* d_c_dQi[j])+ d_v_dQi[j]); 
	       d_fvl_m[3][j] = d_fvl_m[0][j] * ((-(gamma-1.)*u_bar*u_bar+ sign_pm*2.*(gamma-1.)*u_bar*c+2.*c*c)/(gamma*gamma-1.)+.5*(u*u+v*v)) + fvl_m[0] * ((-2.*(gamma-1.)*u_bar*d_ubar_dQi[j]+sign_pm*2.*(gamma-1.)*(d_ubar_dQi[j]*c+u_bar*d_c_dQi[j])+ 4.*c*d_c_dQi[j])/(gamma*gamma-1.)+ u*d_u_dQi[j] + v*d_v_dQi[j]); 
	  }

     }

     /* clean - up*/  
     free(d_rho_dQi);
     free(d_u_dQi);
     free(d_v_dQi);
     free(d_e_dQi);
     free(d_P_dQi);
     free(d_c_dQi);
     free(d_ubar_dQi);

     //completed successfully!
     return 0;
}
