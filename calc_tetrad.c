#include "panhead.h"

void calc_e(double e_lf[4][4], double w_lf[4][4], double g_up[4][4],
	    double g_dn[4][4], double part_p[], double part_v[])
{
  double el,k0,k1,k2,k3,N1,N2,N3;
  int i,j,k;

  el = part_p[3]/part_p[0];
  k1 = part_v[1]/(part_v[0]+part_v[3]*el);
  k2 = (g_up[0][0]+2.*g_up[0][3]*el+g_up[3][3]*el*el)/g_up[1][1];
  k3 = 1.*(part_v[0]+k1*k2*part_v[1]+el*part_v[3]);
  N1 = sqrt(g_dn[0][0]*el*el-2.*g_dn[0][3]*el+g_dn[3][3]);
  N2 = sqrt(g_up[0][0]*k1*k1+2.*g_up[0][3]*k1*k1*el+
	    g_up[1][1]+g_up[3][3]*k1*k1*el*el);
  N3 = sqrt(pow(part_v[2],2.)*(g_up[0][0]+2.*g_up[0][3]*el+
				  g_up[1][1]*k1*k1*k2*k2+g_up[3][3]*el*el)+
	    g_up[2][2]*k3*k3);

  for (j=0;j<=3;j++) e_lf[0][j]=part_v[j];
  e_lf[1][0] = (-g_up[0][0]*k1-el*g_up[0][3]*k1)/N2;
  e_lf[1][1] = g_up[1][1]/N2;
  e_lf[1][2] = 0.;
  e_lf[1][3] = (-g_up[0][3]*k1-g_up[3][3]*el*k1)/N2;
  e_lf[2][0] = (part_v[2]*g_up[0][0]+part_v[2]*g_up[0][3]*el)/N3;
  e_lf[2][1] = g_up[1][1]*k1*k2*part_v[2]/N3;
  e_lf[2][2] = -g_up[2][2]*k3/N3;
  e_lf[2][3] = (part_v[2]*g_up[0][3]+part_v[2]*g_up[3][3]*el)/N3;
  if (N1 > 0) {
    e_lf[3][0] = -el/N1;
    e_lf[3][1] = 0.;
    e_lf[3][2] = 0.;
    e_lf[3][3] = 1./N1;
  } else {
    e_lf[3][0] = 0.;
    e_lf[3][1] = 0.;
    e_lf[3][2] = 0.;
    e_lf[3][3] = 1.;
  }  
  for (i=0;i<=3;i++) {
    for (j=0;j<=3;j++) {
      w_lf[i][j]=0;
      for (k=0;k<=3;k++) {
	w_lf[i][j]+=g_dn[k][i]*e_lf[j][k];
	//w_lf[i][j]+=g_dn[k][i]*g_up[j][k];
      }
    }
  }
}

//calc_e2 creates a right-handed coordinate system with e_x=e_r, e_y=e_phi, and e_z=-e_theta
//when (x,y,z)=(r,0,0) 
void calc_e2(double e_lf[4][4], double w_lf[4][4], double g_up[4][4],
	     double g_dn[4][4], double part_p[], double part_v[])
     
{
  double el,k0,k1,k2,k3,N1,N2,N3;
  int i,j,k;

  el = part_p[3]/part_p[0];
  k1 = part_v[1]/(part_v[0]+part_v[3]*el);
  k2 = (g_up[0][0]+2.*g_up[0][3]*el+g_up[3][3]*el*el)/g_up[1][1];
  k3 = 1.*(part_v[0]+k1*k2*part_v[1]+el*part_v[3]);
  N1 = sqrt(g_dn[0][0]*el*el-2.*g_dn[0][3]*el+g_dn[3][3]);
  N2 = sqrt(g_up[0][0]*k1*k1+2.*g_up[0][3]*k1*k1*el+
	    g_up[1][1]+g_up[3][3]*k1*k1*el*el);
  N3 = sqrt(pow(part_v[2],2.)*(g_up[0][0]+2.*g_up[0][3]*el+
				  g_up[1][1]*k1*k1*k2*k2+g_up[3][3]*el*el)+
	    g_up[2][2]*k3*k3);

  for (j=0;j<=3;j++) e_lf[0][j]=part_v[j];
  e_lf[1][0] = (-g_up[0][0]*k1-el*g_up[0][3]*k1)/N2;
  e_lf[1][1] = g_up[1][1]/N2;
  e_lf[1][2] = 0.;
  e_lf[1][3] = (-g_up[0][3]*k1-g_up[3][3]*el*k1)/N2;
  e_lf[3][0] = (part_v[2]*g_up[0][0]+part_v[2]*g_up[0][3]*el)/N3;
  e_lf[3][1] = g_up[1][1]*k1*k2*part_v[2]/N3;
  e_lf[3][2] = -g_up[2][2]*k3/N3;
  e_lf[3][3] = (part_v[2]*g_up[0][3]+part_v[2]*g_up[3][3]*el)/N3;
  if (N1 > 0) {
    e_lf[2][0] = -el/N1;
    e_lf[2][1] = 0.;
    e_lf[2][2] = 0.;
    e_lf[2][3] = 1./N1;
  } else {
    e_lf[2][0] = 0.;
    e_lf[2][1] = 0.;
    e_lf[2][2] = 0.;
    e_lf[2][3] = 1.;
  }  
  for (i=0;i<=3;i++) {
    for (j=0;j<=3;j++) {
      w_lf[i][j]=0;
      for (k=0;k<=3;k++) {
	w_lf[i][j]+=g_dn[k][i]*e_lf[j][k];
	//w_lf[i][j]+=g_dn[k][i]*g_up[j][k];
      }
    }
  }
}

//calc_e3 creates a right-handed coordinate system in the ZAMO frame

void calc_e3(double e_lf[4][4], double w_lf[4][4], double g_up[4][4],
	     double g_dn[4][4], double part_p[], double part_v[], double y[])
{
  double r,r2,a2,cth,Sig,Delta,alpha,omgbar2;
  int i,j,k;

  r = y[1];
  r2 = r*r;
  cth = cos(y[2]);
  a2 = aa*aa;
  Sig = r2+a2*cth*cth;
  Delta = r2-2*M*r+a2;
  alpha = sqrt(Sig*Delta/(Sig*Delta+2.*M*r*(a2+r2)));
  omgbar2 = Delta/(alpha*alpha)*(1.-cth*cth);
  part_p[0] = -alpha;
  part_p[1] = 0;
  part_p[2] = 0;
  part_p[3] = 0;
  part_v[0] = g_up[0][0]*part_p[0]+g_up[0][3]*part_p[3];
  part_v[1] = g_up[1][1]*part_p[1];
  part_v[2] = g_up[2][2]*part_p[2];
  part_v[3] = g_up[0][3]*part_p[0]+g_up[3][3]*part_p[3];

  for (j=0;j<=3;j++) e_lf[0][j]=part_v[j];
  for (i=1;i<=3;i++) {
    for (j=0;j<=3;j++) e_lf[i][j]=0;
  }
  e_lf[1][1] = sqrt(Delta/Sig);
  e_lf[2][2] = sqrt(1./Sig);
  e_lf[3][3] = sqrt(1./omgbar2);
  for (i=0;i<=3;i++) {
    for (j=0;j<=3;j++) {
      w_lf[i][j]=0;
      for (k=0;k<=3;k++) {
	w_lf[i][j]+=g_dn[k][i]*e_lf[j][k];
      }
    }
  }
}

/*calc_e4 creates a right-handed coordinate system on the surface of the 
photosphere with e_t aligned with the fluid velocity, e_z aligned with 
the normal to the photosphere surface, e_x in the radial direction, and 
e_y in the (+/-) azimuthal direction, depending on the direction of e_z, 
to ensure a right-handed coordinate system.
G_fact is returned as the proper area of the surface patch defined by
dx_r and dx_p.
*/

void calc_e4(double e_lf[4][4], double w_lf[4][4], double g_up[4][4],
	     double g_dn[4][4], double part_p[], double part_v[],
	     double dx_r[], double dx_p[], int ibottom, double *G_fact)
     
{
  double e_t[4],e_x[4],e_y[4],e_z[4],
    e_tlf[4],e_xlf[4],e_ylf[4],e_zlf[4],
    e_x3[3],e_y3[3],e_z3[3],dotx,doty,dd,dr,dph,
    *adata,*bdata,**adata1,*bdata1,*xdata1;
  int i,j,k,*indx;

  indx = (int *)malloc(5*sizeof(double));
  bdata = (double *)malloc(4*sizeof(double));
  adata = (double *)malloc(4*4*sizeof(double));
  xdata1 = (double *)malloc(5*sizeof(double));
  bdata1 = (double *)malloc(5*sizeof(double));
  adata1 = calloc(5,sizeof(double *));
  for (i=0;i<5;i++) adata1[i]=calloc(5,sizeof(double));

  dr = sqrt(dx_r[1]*dx_r[1]+dx_r[2]*dx_r[2]);
  dph = sqrt(dx_p[2]*dx_p[2]+dx_p[3]*dx_p[3]);
  for (i=0;i<=3;i++) e_t[i]=part_v[i];
  dotx = dot_g4(g_dn,dx_r,e_t);
  doty = dot_g4(g_dn,dx_p,e_t);
  for (i=0;i<=3;i++) {
    e_x[i]=dx_r[i]+e_t[i]*dotx;
    e_y[i]=dx_p[i]+e_t[i]*doty;
  }
  //  printf("%12.5e %12.5e %12.5e\n",
  // dot_g4(g_dn,e_x,e_x),dot_g4(g_dn,e_y,e_y),dot_g4(g_dn,e_t,e_t));
  //printf("%12.5e %12.5e %12.5e %12.5e\n",
  //	 e_t[0],e_t[1],e_t[2],e_t[3]);
  //printf("%12.5e %12.5e %12.5e %12.5e\n",
  //	 e_x[0],e_x[1],e_x[2],e_x[3]);
  //printf("%12.5e %12.5e %12.5e %12.5e\n",
  //	 e_y[0],e_y[1],e_y[2],e_y[3]);

  //transform into local gas frame
  calc_e2(e_lf,w_lf,g_up,g_dn,part_p,part_v);
  /*
  printf("\n");
  printf("%12.5e %12.5e %12.5e %12.5e\n",
	 e_lf[0][0],e_lf[0][1],e_lf[0][2],e_lf[0][3]);
  printf("%12.5e %12.5e %12.5e %12.5e\n",
	 e_lf[1][0],e_lf[1][1],e_lf[1][2],e_lf[1][3]);
  printf("%12.5e %12.5e %12.5e %12.5e\n",
	 e_lf[2][0],e_lf[2][1],e_lf[2][2],e_lf[2][3]);
  printf("%12.5e %12.5e %12.5e %12.5e\n",
	 e_lf[3][0],e_lf[3][1],e_lf[3][2],e_lf[3][3]);
  */
  for (i=0;i<=3;i++) {
    for (j=0;j<=3;j++) {
      adata1[i+1][j+1]=e_lf[j][i];
    }
  }
  
  ludcmp_js(adata1,4,indx,&dd);
  for (j=0;j<=3;j++) bdata1[j+1]=e_t[j];
  lubksb_js(adata1,4,indx,bdata1);
  for (j=0;j<=3;j++) e_tlf[j]=bdata1[j+1];
  for (j=0;j<=3;j++) bdata1[j+1]=e_x[j];
  lubksb_js(adata1,4,indx,bdata1);
  for (j=0;j<=3;j++) e_xlf[j]=bdata1[j+1];
  for (j=0;j<=3;j++) bdata1[j+1]=e_y[j];
  lubksb_js(adata1,4,indx,bdata1);
  for (j=0;j<=3;j++) e_ylf[j]=bdata1[j+1];
  for (j=0;j<=3;j++) e_zlf[j]=0;
  
  for (i=0;i<=2;i++) {
    e_x3[i]=e_xlf[i+1];
    e_y3[i]=e_ylf[i+1];
  }
  cross(e_x3,e_y3,e_z3);
  //printf("ez3 %12.5e\n",calc_mag(e_z3));
  *G_fact = calc_mag(e_z3);
  normalize(e_x3);
  normalize(e_y3);
  normalize(e_z3);
  cross(e_z3,e_x3,e_y3);
  //top side of disk: e_z = e_r x e_ph
  if (ibottom == 0) {
    for (i=0;i<=2;i++) {
      e_zlf[i+1]=e_z3[i];
      e_ylf[i+1]=e_y3[i];
      e_xlf[i+1]=e_x3[i];
    }
  }
  //bottom side of disk: e_z = -e_r x e_ph; e_x = e_r; e_y = -e_ph
  if (ibottom == 1) {
    for (i=0;i<=2;i++) {
      e_zlf[i+1]=-e_z3[i];
      e_ylf[i+1]=-e_y3[i];
      e_xlf[i+1]=e_x3[i];
    }
  }
  /*
  printf("e_tlf: %12.5e %12.5e %12.5e %12.5e\n",
	 e_tlf[0],e_tlf[1],e_tlf[2],e_tlf[3]);
  printf("e_xlf: %12.5e %12.5e %12.5e %12.5e\n",
	 e_xlf[0],e_xlf[1],e_xlf[2],e_xlf[3]);
  printf("e_ylf: %12.5e %12.5e %12.5e %12.5e\n",
	 e_ylf[0],e_ylf[1],e_ylf[2],e_ylf[3]);
  printf("e_zlf: %12.5e %12.5e %12.5e %12.5e\n",
	 e_zlf[0],e_zlf[1],e_zlf[2],e_zlf[3]);
  */
  for (i=0;i<=3;i++) {
    e_t[i]=0;
    e_x[i]=0;
    e_y[i]=0;
    e_z[i]=0;
    for (j=0;j<=3;j++) {
      e_t[i]+=e_lf[j][i]*e_tlf[j];
      e_x[i]+=e_lf[j][i]*e_xlf[j];
      e_y[i]+=e_lf[j][i]*e_ylf[j]; 
      e_z[i]+=e_lf[j][i]*e_zlf[j];
    }
  }
  for (i=0;i<=3;i++) {
    e_lf[0][i]=e_t[i];
    e_lf[1][i]=e_x[i];
    e_lf[2][i]=e_y[i];
    e_lf[3][i]=e_z[i];
  }
  for (i=0;i<=3;i++) {
    for (j=0;j<=3;j++) {
      w_lf[i][j]=0;
      for (k=0;k<=3;k++) {
	w_lf[i][j]+=g_dn[k][i]*e_lf[j][k];
      }
    }
  }
  free(indx);
  free(adata);
  free(bdata);
  free(bdata1);
  free(xdata1);
  for (i=0;i<=4;i++) free(adata1[i]);
  free(adata1);
  /*
  printf("%12.5e %12.5e %12.5e %12.5e\n",
	 e_t[0],e_t[1],e_t[2],e_t[3]);
  printf("%12.5e %12.5e %12.5e %12.5e\n",
	 e_x[0],e_x[1],e_x[2],e_x[3]);
  printf("%12.5e %12.5e %12.5e %12.5e\n",
	 e_y[0],e_y[1],e_y[2],e_y[3]);
  printf("%12.5e %12.5e %12.5e %12.5e\n",
	 e_z[0],e_z[1],e_z[2],e_z[3]);
  printf("%12.2e %12.2e %12.2e %12.2e\n",
	 dot_g4(g_dn,e_t,e_t),dot_g4(g_dn,e_t,e_x),dot_g4(g_dn,e_t,e_y),dot_g4(g_dn,e_t,e_z));
  printf("%12.2e %12.2e %12.2e %12.2e\n",
	 dot_g4(g_dn,e_x,e_t),dot_g4(g_dn,e_x,e_x),dot_g4(g_dn,e_x,e_y),dot_g4(g_dn,e_x,e_z));
  printf("%12.2e %12.2e %12.2e %12.2e\n",
	 dot_g4(g_dn,e_y,e_t),dot_g4(g_dn,e_y,e_x),dot_g4(g_dn,e_y,e_y),dot_g4(g_dn,e_y,e_z));
  printf("%12.2e %12.2e %12.2e %12.2e\n",
	 dot_g4(g_dn,e_z,e_t),dot_g4(g_dn,e_z,e_x),dot_g4(g_dn,e_z,e_y),dot_g4(g_dn,e_z,e_z));
  */
}

/*calc_dV calculates the proper volume element spanned by dx_r, dx_th, 
and dx_p, and returns it as V_fact.
*/

void calc_dV(double e_lf[4][4], double w_lf[4][4], double g_up[4][4],
	     double g_dn[4][4], double part_p[], double part_v[],
	     double dx_r[], double dx_th[], double dx_p[], 
	     double *V_fact, double *G_fact)
     
{
  double e_t[4],e_x[4],e_y[4],e_z[4],
    e_tlf[4],e_xlf[4],e_ylf[4],e_zlf[4],
    e_x3[3],e_y3[3],e_z3[3],n_3[3],dotx,doty,dotz,dd,
    *adata,*bdata,**adata1,*bdata1,*xdata1;
  int i,j,k,*indx;

  indx = (int *)malloc(5*sizeof(double));
  bdata = (double *)malloc(4*sizeof(double));
  adata = (double *)malloc(4*4*sizeof(double));
  xdata1 = (double *)malloc(5*sizeof(double));
  bdata1 = (double *)malloc(5*sizeof(double));
  adata1 = calloc(5,sizeof(double *));
  for (i=0;i<5;i++) adata1[i]=calloc(5,sizeof(double));

  for (i=0;i<=3;i++) e_t[i]=part_v[i];
  dotx = dot_g4(g_dn,dx_r,e_t);
  doty = dot_g4(g_dn,dx_p,e_t);
  dotz = dot_g4(g_dn,dx_th,e_t);
  for (i=0;i<=3;i++) {
    e_x[i]=dx_r[i]+e_t[i]*dotx;
    e_y[i]=dx_p[i]+e_t[i]*doty;
    e_z[i]=dx_th[i]+e_t[i]*dotz;
  }
  //printf("%12.5e %12.5e %12.5e\n",
  // dot_g4(g_dn,e_x,e_x),dot_g4(g_dn,e_y,e_y),dot_g4(g_dn,e_z,e_z));
  //printf("%12.5e %12.5e %12.5e %12.5e\n",
  // e_t[0],e_t[1],e_t[2],e_t[3]);
  //printf("%12.5e %12.5e %12.5e %12.5e\n",
  // e_x[0],e_x[1],e_x[2],e_x[3]);
  //printf("%12.5e %12.5e %12.5e %12.5e\n",
  // e_y[0],e_y[1],e_y[2],e_y[3]);

  //transform into local gas frame
  calc_e2(e_lf,w_lf,g_up,g_dn,part_p,part_v);
  /*
  printf("\n");
  printf("%12.5e %12.5e %12.5e %12.5e\n",
	 e_lf[0][0],e_lf[0][1],e_lf[0][2],e_lf[0][3]);
  printf("%12.5e %12.5e %12.5e %12.5e\n",
	 e_lf[1][0],e_lf[1][1],e_lf[1][2],e_lf[1][3]);
  printf("%12.5e %12.5e %12.5e %12.5e\n",
	 e_lf[2][0],e_lf[2][1],e_lf[2][2],e_lf[2][3]);
  printf("%12.5e %12.5e %12.5e %12.5e\n",
	 e_lf[3][0],e_lf[3][1],e_lf[3][2],e_lf[3][3]);
  */
  for (i=0;i<=3;i++) {
    for (j=0;j<=3;j++) {
      adata1[i+1][j+1]=e_lf[j][i];
    }
  }
  
  ludcmp_js(adata1,4,indx,&dd);
  for (j=0;j<=3;j++) bdata1[j+1]=e_t[j];
  lubksb_js(adata1,4,indx,bdata1);
  for (j=0;j<=3;j++) e_tlf[j]=bdata1[j+1];
  for (j=0;j<=3;j++) bdata1[j+1]=e_x[j];
  lubksb_js(adata1,4,indx,bdata1);
  for (j=0;j<=3;j++) e_xlf[j]=bdata1[j+1];
  for (j=0;j<=3;j++) bdata1[j+1]=e_y[j];
  lubksb_js(adata1,4,indx,bdata1);
  for (j=0;j<=3;j++) e_ylf[j]=bdata1[j+1];
  for (j=0;j<=3;j++) bdata1[j+1]=e_z[j];
  lubksb_js(adata1,4,indx,bdata1);
  for (j=0;j<=3;j++) e_zlf[j]=bdata1[j+1];
  
  for (i=0;i<=2;i++) {
    e_x3[i]=e_xlf[i+1];
    e_y3[i]=e_ylf[i+1];
    e_z3[i]=e_zlf[i+1];
  }
  cross(e_x3,e_y3,n_3);
  *G_fact = calc_mag(n_3);
  *V_fact = fabs(dot(n_3,e_z3));
  for (i=0;i<=2;i++) e_z3[i]=n_3[i];
  normalize(e_x3);
  normalize(e_y3);
  normalize(e_z3);
  cross(e_z3,e_x3,e_y3);
  for (i=0;i<=2;i++) {
    e_zlf[i+1]=e_z3[i];
    e_ylf[i+1]=e_y3[i];
    e_xlf[i+1]=e_x3[i];
  }
  for (i=0;i<=3;i++) {
    e_t[i]=0;
    e_x[i]=0;
    e_y[i]=0;
    e_z[i]=0;
    for (j=0;j<=3;j++) {
      e_t[i]+=e_lf[j][i]*e_tlf[j];
      e_x[i]+=e_lf[j][i]*e_xlf[j];
      e_y[i]+=e_lf[j][i]*e_ylf[j]; 
      e_z[i]+=e_lf[j][i]*e_zlf[j];
    }
  }
  for (i=0;i<=3;i++) {
    e_lf[0][i]=e_t[i];
    e_lf[1][i]=e_x[i];
    e_lf[2][i]=e_y[i];
    e_lf[3][i]=e_z[i];
  }
  for (i=0;i<=3;i++) {
    for (j=0;j<=3;j++) {
      w_lf[i][j]=0;
      for (k=0;k<=3;k++) {
	w_lf[i][j]+=g_dn[k][i]*e_lf[j][k];
      }
    }
  }

  free(indx);
  free(adata);
  free(bdata);
  free(bdata1);
  free(xdata1);
  for (i=0;i<=4;i++) free(adata1[i]);
  free(adata1);
}

void surface_tetrads(double rr[], double tt[], double pp[],
		     double ut_ijk[], double ur_ijk[], double uz_ijk[],
		     double up_ijk[], int diskbody_ik[], 
		     double emtop_ik[], double embot_ik[],
		     double reftop_ik[], double refbot_ik[],
		     double emtop_elf_ik[], double embot_elf_ik[],
		     double reftop_elf_ik[], double refbot_elf_ik[],
		     double rho_ijk[], double T_ijk[], double bb_ijk[],
		     double Gtop_ik[], double Gbot_ik[], 
		     long irstart, double eps_th)
{
  double part_x0[4],part_v[4], part_p[4],g_up[4][4],g_dn[4][4],
    e_lf[4][4],w_lf[4][4],dx_r[4],dx_p[4],weights[12],
    rho0,T_e0,Bmag0,wlo,whi,G_fact;
  long ibottom,ir,iph,jth_lo,jth_hi,kph,i,j;
  for (ibottom=0;ibottom<=1;ibottom++) {
    for (ir=irstart;ir<=Nr;ir++) {
      for (iph=0;iph<=Nph;iph++) {
	part_x0[0]=0;
	part_x0[1]=rr[ir];
	if (ibottom == 0) part_x0[2]=emtop_ik[indexr(ir,iph)]-eps_th;
	if (ibottom == 1) part_x0[2]=embot_ik[indexr(ir,iph)]+eps_th;
	part_x0[3]=pp[iph];
	calc_g(g_dn,g_up,part_x0);
	lookup_data(part_x0,rr,tt,pp,g_dn,
		    rho_ijk,T_ijk,bb_ijk,ut_ijk,ur_ijk,uz_ijk,up_ijk,
		    weights,&rho0,&T_e0,&Bmag0,part_v);
	part_p[0] = g_dn[0][0]*part_v[0]+g_dn[0][3]*part_v[3];
	part_p[1] = g_dn[1][1]*part_v[1];
	part_p[2] = g_dn[2][2]*part_v[2];
	part_p[3] = g_dn[3][0]*part_v[0]+g_dn[3][3]*part_v[3];
	
	dx_r[0] = 0;
	dx_r[1] = 0.5*(rr[ir+1]-rr[ir-1]);
	dx_r[2] = 0;
	dx_r[3] = 0;
	if ((diskbody_ik[indexr(ir,iph)] > 0)&&
	    (diskbody_ik[indexr(ir+1,iph)] > 0)&&
	    (diskbody_ik[indexr(ir-1,iph)] > 0)) {
	  if (ibottom == 0) 
	    dx_r[2] = 0.5*(emtop_ik[indexr(ir+1,iph)]-
			   emtop_ik[indexr(ir-1,iph)]);
	  if (ibottom == 1) 
	    dx_r[2] = 0.5*(embot_ik[indexr(ir+1,iph)]-
			   embot_ik[indexr(ir-1,iph)]);
	}
	dx_p[0] = 0;
	dx_p[1] = 0;
	dx_p[2] = 0;
	dx_p[3] = (pp[1]-pp[0]);
	if ((diskbody_ik[indexr(ir,iph-1)] > 0)&&
	    (diskbody_ik[indexr(ir,iph)] > 0)&&
	    (diskbody_ik[indexr(ir,iph+1)] > 0)) {
	  if (ibottom == 0)
	    dx_p[2] = 0.5*(emtop_ik[indexr(ir,iph+1)]-
			   emtop_ik[indexr(ir,iph-1)]);
	  if (ibottom == 1)
	    dx_p[2] = 0.5*(embot_ik[indexr(ir,iph+1)]
			   -embot_ik[indexr(ir,iph-1)]);
	}
	calc_e4(e_lf,w_lf,g_up,g_dn,part_p,part_v,dx_r,dx_p,ibottom,&G_fact);
	if (ibottom == 0) {
	  for (i=0;i<=3;i++) {
	    for (j=0;j<=3;j++) {
	      emtop_elf_ik[indexelf(ir,iph,i,j)]=e_lf[i][j];
	    }
	  }
	  Gtop_ik[indexr(ir,iph)]=G_fact;
	}
	if (ibottom == 1) {
	  for (i=0;i<=3;i++) {
	    for (j=0;j<=3;j++) {
	      embot_elf_ik[indexelf(ir,iph,i,j)]=e_lf[i][j];
	    }
	  }
	  Gbot_ik[indexr(ir,iph)]=G_fact;
	}
	//printf("%d %d %d %10.4e %10.4e %10.4e %10.4e\n",
	//     ibottom,ir,iph,part_v[0],part_v[1],part_v[2],part_v[3]);
	//printf("%10.4e %10.4e %10.4e %10.4e\n",
	//     e_lf[0][0],e_lf[0][1],e_lf[0][2],e_lf[0][3]);
	dx_r[0] = 0;
	dx_r[1] = 0.5*(rr[ir+1]-rr[ir-1]);
	dx_r[2] = 0;
	dx_r[3] = 0;
	if ((diskbody_ik[indexr(ir,iph)] > 1)&&
	    (diskbody_ik[indexr(ir+1,iph)] > 1)&&
	    (diskbody_ik[indexr(ir-1,iph)] > 1)) {
	  if (ibottom == 0) 
	    dx_r[2] = 0.5*(reftop_ik[indexr(ir+1,iph)]-
			   reftop_ik[indexr(ir-1,iph)]);
	  if (ibottom == 1) 
	    dx_r[2] = 0.5*(refbot_ik[indexr(ir+1,iph)]-
			   refbot_ik[indexr(ir-1,iph)]);
	}
	dx_p[0] = 0;
	dx_p[1] = 0;
	dx_p[2] = 0;
	dx_p[3] = (pp[1]-pp[0]);
	if ((diskbody_ik[indexr(ir,iph-1)] > 1)&&
	    (diskbody_ik[indexr(ir,iph)] > 1)&&
	    (diskbody_ik[indexr(ir,iph+1)] > 1)) {
	  if (ibottom == 0)
	    dx_p[2] = 0.5*(reftop_ik[indexr(ir,iph+1)]-
			   reftop_ik[indexr(ir,iph-1)]);
	  if (ibottom == 1)
	    dx_p[2] = 0.5*(refbot_ik[indexr(ir,iph+1)]
			   -refbot_ik[indexr(ir,iph-1)]);
	}
	calc_e4(e_lf,w_lf,g_up,g_dn,part_p,part_v,dx_r,dx_p,ibottom,&G_fact);
	if (ibottom == 0) {
	  for (i=0;i<=3;i++) {
	    for (j=0;j<=3;j++) {
	      reftop_elf_ik[indexelf(ir,iph,i,j)]=e_lf[i][j];
	    }
	  }
	}
	if (ibottom == 1) {
	  for (i=0;i<=3;i++) {
	    for (j=0;j<=3;j++) {
	      refbot_elf_ik[indexelf(ir,iph,i,j)]=e_lf[i][j];
	    }
	  }
	}
      }
    }
  }  
}

