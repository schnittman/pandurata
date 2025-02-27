#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include "mpi.h"

#define M 1.0
#define Mstar 3.33
//#define Mdot_17 8.3
#define L_Edd 0.01
//scale_height = 1.0
//surface and Temp smoothed
#define del_eta 0.0
#define nt_alpha 0.1
#define aa 0.0
#define tau_es 100.0
#define T_cor 100.0 //corona temp in keV
#define scale_height 0.2
//1 for line emission, 2 for N-T thermal blackbody, 3 for R^{-3/4} thermal BB
#define atm_model 4
#define Nclumps 10000
#define overdensity 10.0
#define OPT_THIN 0.0
#define TWO_SIDED 1.0
#define em_model 2.0
#define RUN_ID 1251.0

#define N 10
#define Ne 100
#define Nr 179
#define Nth 159
#define Nph 255
#define Ni 200
#define Nt 20
#define Nph_obs 39
#define Nth_obs 40
#define Ne_obs 10
#define Ns_obs 5
#define N_symm 1

#define PI 3.14159265358979323846
#define kB 1.38e-16
#define kB_ev 8.6173e-5
#define ak 7.566e-25
#define k_es 0.4
#define sigma 5.67e-5
#define sigma_T 6.65e-25
#define cc 3.0e10
#define Gn 6.67e-8
#define mp 1.67e-24
#define me 9.11e-28
#define qe 4.803e-10

/******************INDEXING MACROS*******************************/
#define indexi(a,b,c) ((Nth_obs+1)*(Ni+1)*(c)+(Nth_obs+1)*(b)+(a))
#define indexspci(a,b,c,d) ((Nth_obs+1)*(Ni+1)*(Ni+1)*(d)+(Nth_obs+1)*(Ni+1)*(c)+(Nth_obs+1)*(b)+(a))
#define indexphi(a,b,c,d) ((Nth_obs+1)*(Nph_obs+1)*(Ni+1)*(d)+(Nth_obs+1)*(Nph_obs+1)*(c)+(Nth_obs+1)*(b)+(a))
#define indexrth(a,b,c) ((Nr+1)*(Nth_obs+1)*(c)+(Nr+1)*(b)+(a))
#define index3(a,b,c) ((Nth_obs+1)*(Ne+1)*(c)+(Nth_obs+1)*(b)+(a))
#define indexph(a,b,c) ((Nth_obs+1)*(Nph_obs+1)*(c)+(Nth_obs+1)*(b)+(a))
#define indexijk(a,b,c) ((Nph+1)*(Nth+1)*(a)+(Nph+1)*(b)+(c))
#define indexthijk(a,b,c) ((Nph+1)*(Nth+2)*(a)+(Nph+1)*(b)+(c))
#define indexelf(a,b,c,d) ((Nph+1)*(4)*(4)*(a)+(4)*(4)*(b)+(4)*(c)+(d))
#define index2(a,b) ((Nth_obs+1)*(b)+(a))
#define indexr(a,b) ((Nr+1)*(int)(fmod(b+Nph+1,Nph+1))+(a))
#define indexre(a,b) ((Nr+1)*(b)+(a))
#define indexrpe(a,b,c) ((Nr+1)*(Nph+1)*(c)+(Nr+1)*(b)+(a))
#define indexc(a,b) ((Nclumps+1)*(b)+(a))

/****************FUNCTION PROTOTYPES*****************************/
double F_sync(double x);
double F_cycl(double x, double T_e);
double lookup_Pnorm(double Pnorm[], double T_[], double T_e, int N_T);
void calc_Pnorm_T(double Pnorm[], double T[], int N_T);
void m_mult3(double A[3][3], double B[3][3], double C[3][3]);
double ran2(long *idum);
  

/****************VECTOR_MATH.C***********************************/
double calc_mag(double x[]);
void cross(double a[], double b[], double c[]);
void normalize(double x[]);
double dot(double a[], double b[]);
void cart_to_spher(double x[], double r[]);
void spher_to_cart(double r[], double x[]);
void cartv_to_spherv(double x[], double v[], double p[]);
void spherv_to_cartv(double r[], double p[], double v[]);

/****************TENSOR_MATH.C***********************************/
double dot_g4(double g_ab[4][4],double a[],double b[]);
void norml_tl(double g_ab[4][4], double a[]);
void norml_sl(double g_ab[4][4], double a[]);
void boost(double beta, double n_[], double p_[]);
void ludcmp_js(double **a, int n, int *indx, double *d);
void lubksb_js(double **a, int n, int *indx, double b[]);

/****************CALC_G.C****************************************/
void calc_g(double g_dn[4][4], double g_up[4][4], double y[]);

/****************CALC_TETRAD.C***********************************/
void calc_e(double e_lf[4][4], double w_lf[4][4], double g_up[4][4],
	    double g_dn[4][4], double part_p[], double part_v[]);
void calc_e2(double e_lf[4][4], double w_lf[4][4], double g_up[4][4],
	     double g_dn[4][4], double part_p[], double part_v[]);
void calc_e3(double e_lf[4][4], double w_lf[4][4], double g_up[4][4],
	     double g_dn[4][4], double part_p[], double part_v[], double y[]);
void calc_e4(double e_lf[4][4], double w_lf[4][4], double g_up[4][4],
	     double g_dn[4][4], double part_p[], double part_v[],
	     double dx_r[], double dx_p[], int ibottom, double *G_fact);
void calc_dV(double e_lf[4][4], double w_lf[4][4], double g_up[4][4],
	     double g_dn[4][4], double part_p[], double part_v[],
	     double dx_r[], double dx_th[], double dx_p[], 
	     double *V_fact, double *G_fact);
void surface_tetrads(double rr[], double tt[], double pp[],
		     double ut_ijk[], double ur_ijk[], double uz_ijk[],
		     double up_ijk[], int diskbody_ik[], 
		     double emtop_ik[], double embot_ik[],
		     double reftop_ik[], double refbot_ik[],
		     double emtop_elf_ik[], double embot_elf_ik[],
		     double reftop_elf_ik[], double refbot_elf_ik[],
		     double rho_ijk[], double T_ijk[], double bb_ijk[],
		     double Gtop_ik[], double Gbot_ik[], 
		     long irstart, double eps_th);
  
/****************CHANDRA.C**************************************/
void chandra_limit(double mu, double *deg, double *darken);
void diffuse_reflection(double mu0, double ph0, double mu, double ph,
			double *I_l, double *I_r, double *U_);

void nt_spectrum(double Risco, double rr[], double drr[], 
		 double nu[], double Inur[], double Ts_r[], double qnur[]);
void get_harm3d_data(double rr[], double tt[], double pp[],
		     double rho_ijk[], double T_ijk[], double bb_ijk[], double tau_ijk[],
		     double ut_ijk[], double ur_ijk[], double uz_ijk[], 
		     double up_ijk[], int diskbody_ik[], double sigtau_ik[], 
		     double Tdisk_ik[], double emtop_ik[], double embot_ik[], 
		     double reftop_ik[], double refbot_ik[]);
void lookup_data(double part_x0[], 
		 double rr[], double tt[], double pp[], double g_dn[4][4],
		 double rho_ijk[], double T_ijk[], double bb_ijk[], double ut_ijk[],
		 double ur_ijk[], double uz_ijk[], double up_ijk[],
		 double weights[], double *rho0, double *T_e0, double *Bmag0, double part_v[]);
void calc_scat_angles(double deg, double *mu, double *psi);
void cashkarp(double dt, double y0[], double yn[], double del[]);
void accel(double y[], double dy[]);
void advance(double dt, double y0[], double yn[]);

/****************TIME_KEEPER.C**********************************/
void start_time(void);
double wrt_time(void);
