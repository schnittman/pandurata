#include "panhead.h"
#include <unistd.h>

#define erro 1e-8
#define emin 1e-14

double F_sync(double x)
{
  double Ftot,Flo,Fhi,gam_13;
  gam_13 = 2.67894;
  Flo = 4*PI/sqrt(3.)/gam_13*pow(x/2.,1./3.)*
    (1.-gam_13/2.*pow(x/2.,2./3.)+3./4.*(x*x/4.));
  Fhi = sqrt(PI/2.)*exp(-x)*sqrt(x)*(1.+55./72./x);
  Ftot = Flo*(exp(-x*x))+Fhi*(1.-exp(-x*x));
  Ftot = Ftot/1.138; //fix normalization
  return Ftot;
}

//T_e corona temp in keV
double F_cycl(double x, double T_e)
{
  double Ftot;
  Ftot = sqrt(1./x-1.)*exp(-511./T_e*(1./x-1.));
  return Ftot;
}

double lookup_Pnorm(double Pnorm[], double T_[], double T_e, int N_T)
{
  double wlo, whi, PN_T;
  int ilo, ihi;
  ilo=(int)floor((N_T)*log(T_e/T_[0])/log(T_[N_T]/T_[0]));
  if (ilo < 0) PN_T = Pnorm[0];
  if (ilo >= N_T) PN_T = Pnorm[N_T];
  if ((ilo >= 0)&&(ilo < N_T)) {
    ihi = ilo+1;
    wlo = (T_[ihi]-T_e)/(T_[ihi]-T_[ilo]);
    whi = (T_e-T_[ilo])/(T_[ihi]-T_[ilo]);
    PN_T = Pnorm[ilo]*wlo+Pnorm[ihi]*whi;
  }
  //printf("lookup %g %g\n",T_e,PN_T);
  return PN_T;
}

//tabulate the normalization of cyclotron emission as a function of 
//electron temperature T (measured in keV)
void calc_Pnorm_T(double Pnorm[], double T[], int N_T)
{
  int i, i_t, Nx=1001;
  double x_lo, x_hi, x0,x1,dlogx, dx, x, Px, Ptot;
  x_lo = 1e-3;
  x_hi = 1;
  dlogx = log10(x_hi/x_lo)/(Nx-1.);
  for (i_t=0;i_t<N_T;i_t++) {
    x0 = x_lo;
    Ptot = 0;
    for (i=1;i<Nx;i++) {
      x1 = x_lo*pow(10.,i*dlogx);
      dx = x1-x0;
      x = 0.5*(x1+x0);
      Px = sqrt(1./x-1.)*exp(-511./T[i_t]*(1./x-1.));
      Ptot+=Px*dx;
      //      printf("%10.5e %10.5e %10.5e\n",x,dx,Px);
      x0 = x1;
    }
    if (Ptot < 1e-100) Ptot = 1e-10;
    Pnorm[i_t]=Ptot;
    //printf("%10.5e %10.5e\n",T[i_t],Pnorm[i_t]);
  }
}

void m_mult3(double A[3][3], double B[3][3], double C[3][3])
{
  int i,j,k;
  for (i=0;i<=2;i++) {
    for (j=0;j<=2;j++) {
      C[i][j]=0;
      for (k=0;k<=2;k++) {
	C[i][j]+=A[i][k]*B[k][j];
      }
    }
  }
}


int main(int argc, char* argv[])
{
  double t,r,dth,dph,f,dt,cth,sth,cth2,sth2,phi,ph0,ph2,mu,mu0,err,x,
    I_ll,I_rr,U_,f_new,f_old,dcth,Ncth,th0,wc,nuc,Bmag,gam2,P0,dl,
    y[8],y1[8],y2[8],yn[8],del[8],y_ck[8],del_ck[8],e_x[3],e_y[3],
    x_[3],p_0[3],p_x[3],p_y[3],p_hat[3],r_hat[3],n_hat[3],f_hat[3],fp_hat[3],z_hat[3],
    part_x0[4],part_x[4],part_p[4],part_v[4],part_v_hat[4],f0_[4],f_i[4],f_[4],
    ph_p[4],ph_v[4],ph_v_hat[4],ph_v_p[4],ph_p_p[4],v_[4],p_[4],po_[4],k_[4],
    f_v_hat[4],dx_r[4],dx_th[4],dx_p[4],
    g_up[4][4],g_dn[4][4],lambda_hat[4][4],lam2_hat[4][4],
    g_up_ph[4][4],g_dn_ph[4][4],e_lf[4][4],w_lf[4][4],e_lfs[4][4],w_lfs[4][4],
    e_lf2[4][4],w_lf2[4][4],
    e_perp[3],e_parl[3],e_perp_f[3],e_parl_f[3],e_x_hat[3],e_y_hat[3],
    e_z_hat[3],n_p_hat[3],
    nu0[Ne+1],nu[Ne+1],dnu0[Ne+1],dnu[Ne+1],rr[Nr+1],drr[Nr+1],rrb[Nr+2],
    ntA[Nr+1],ntB[Nr+1],ntC[Nr+1],ntD[Nr+1],ntE[Nr+1],ntF[Nr+1],ntG[Nr+1],
    ntT[Nr+1],nui0[Ne_obs+1],dnui0[Ne_obs+1],
    tt[Nth+1],pp[Nph+1],dtt[Nth+1],dpp[Nph+1],weights[12],wght,
    Risco,Z1,Z2,R_g,tau_tot,tau_tt,l_tot,za,zb,zc,s1,s2,zq,mu_p,beta,gamma,cpsi,spsi,
    T_e,T_e0,T_d,T_d0,rho0,rho,n_e,rdsh,rdsh3,f_hard,flux,flux2,
    kap10,kap20,kap1i,kap2i,normf,cth0,deg,psi,Rout_flux,
    e_min,e_max,phase,tp,E_i,E_f,E_abs,E_em,Rin,Rring,Rout,Redge,dR,A_fact,A_fact2,
    L_fact,B_fact,G_fact,T_fact,D_fact,V_fact,V_tot,h_fact,
    I_p,U_p,Q_p,psip,degp,Omg_d,Eph,f_Fe,N_Fe,sig_edge,
    Ecirc,Lcirc,Eisco,Lisco,rdot,phidot,pow_i,pow_f,rdsh0,
    pro,Rshell,dtau,dtau_es,kapp_es,r_es,rho_bar,R_atm,dR_atm,
    z_screen,delmag,y2mag,E0,pdv_em,Sig,Delta,Carter,lambda,
    alpha,omgb,omg,omgk,alpha2,omgb2,
    Carter0,a2,r2,r3,Rhor,wlo,whi,Nescape,Ndisk,Nscat,Ncapt,th_i,fov,dd,
    Tmin,Tmax,Tobs,dTbin,mfp,dt_min,Iptot,Ipprod,eps_th,fcseed,
    Afact_[Ne+1],Afact2_[Ne+1],Bnu_[Ne+1],B_nu[Ne+1],j_nu[Ne+1],
    alpha_nu[Ne+1],tau_nu[Ne+1],E0_[Ne+1],Iph_[Ne+1],T_tab[Ne+1],Pnorm[Ne+1],
    *rho_ijk,*T_ijk,*bb_ijk,*tau_ijk,*ut_ijk,*ur_ijk,*uz_ijk,*up_ijk,
    *sigtau_ik,*Tdisk_ik,*emtop_ik,*embot_ik,*reftop_ik,*refbot_ik,
    *emtop_elf_ik,*embot_elf_ik,*reftop_elf_ik,*refbot_elf_ik,
    *Gtop_ik,*Gbot_ik,*corpow_ijk,
    *y_print,*spec,*Iobs,*Ispec,*Qspec,*Uspec,*Ispec_s,*Qspec_s,*Uspec_s,
    *Ispecr,*Rspecr,*Qspecr,*Uspecr,*Cspecr,*Lspec,*Inur,*Inur_NT,*qnur,*Rspec,
    *movie,*y_pass,*l_pass,*I_r,*I_rph,*I_pass,*R_pass,*Ispecp,
    *L_factr,*G_factr,*B_factr,*T_factr,
    *adata,*bdata,**adata1,*bdata1,*xdata1,*image,*imagex,*imagey,
    *spcimage,*spcimagex,*spcimagey,*phimage,*phimagex,*phimagey,
    *x_clumps,*r_clumps;
  double y0[8]={0,10,0,0,1,0,0,0},nu_cut=500.,dnu_cut=100.;
  long it,ip,ith,iph,jth,jph,je,jt,i,j,k,iv,ix,iy,ir,jr,Fe_phot,elo,ehi,
    view_jph,jquad,jphq,iv_Fe,
    file_ik,file_ih,file_iu,file_id,Nscat_tot[100],
    acc,steps,raydex,moddex,Nbins,dscat,iscat,isort,bigprint,imageprint,
    spec_model,rr_model,c_frame,just_scattered,irstart,irstop,irstep,
    iphstart,iphstop,iphstep;
  int signum,in_clump,ir_old,ir_new,jph_old,jph_new,ibottom,isbottom,
    old_zone,new_zone,*indx,*diskbody_ik;
  double big,temp;
  FILE *outfile,*outfile2;
  char fname1[100],fname2[100],fname3[100],fname4[100],fname5[100],fname6[100],
    fname7[100],fname8[100],fname9[100];

  int myid,numprocs,tag,source;
  /*
  MPI_Status status;
  int comm = MPI_COMM_WORLD;
  MPI_Datatype doubletype=MPI_DOUBLE;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank (MPI_COMM_WORLD,&myid);
  MPI_Comm_size (MPI_COMM_WORLD,&numprocs);
  tag=1;
  */

  start_time();
  bigprint = 0;
  imageprint = 1;
  myid = 0;
  numprocs = 1;
  indx = (int *)malloc(5*sizeof(double));
  bdata = (double *)malloc(4*sizeof(double));
  adata = (double *)malloc(4*4*sizeof(double));
  xdata1 = (double *)malloc(5*sizeof(double));
  bdata1 = (double *)malloc(5*sizeof(double));
  adata1 = calloc(5,sizeof(double *));
  for (i=0;i<5;i++) adata1[i]=calloc(5,sizeof(double));
  Iobs = (double *)malloc((Nth_obs+1)*(Nt+1)*(Nph_obs+1)*sizeof(double));
  Ispec = (double *)malloc((Nth_obs+1)*(Ne+1)*sizeof(double));
  Qspec = (double *)malloc((Nth_obs+1)*(Ne+1)*sizeof(double));
  Uspec = (double *)malloc((Nth_obs+1)*(Ne+1)*sizeof(double));
  Ispec_s = (double *)malloc((Nth_obs+1)*(Ne+1)*6*sizeof(double));
  Qspec_s = (double *)malloc((Nth_obs+1)*(Ne+1)*6*sizeof(double));
  Uspec_s = (double *)malloc((Nth_obs+1)*(Ne+1)*6*sizeof(double));
  Ispecr = (double *)malloc((Nr+1)*(Nth_obs+1)*(Ne+1)*sizeof(double));
  Rspecr = (double *)malloc((Nr+1)*(Nth_obs+1)*(Ne+1)*sizeof(double));
  Qspecr = (double *)malloc((Nr+1)*(Nth_obs+1)*(Ne+1)*sizeof(double));
  Uspecr = (double *)malloc((Nr+1)*(Nth_obs+1)*(Ne+1)*sizeof(double));
  Cspecr = (double *)malloc((Nr+1)*(Ne+1)*sizeof(double));
  Ispecp = (double *)malloc((Nth_obs+1)*(Nph_obs+1)*(Ne+1)*sizeof(double));
  Lspec = (double *)malloc((Nth_obs+1)*(Ne+1)*sizeof(double));
  Rspec = (double *)malloc((Nr+1)*(Ne+1)*sizeof(double));
  Inur = (double *)malloc((Nr+1)*(Ne+1)*sizeof(double));
  Inur_NT = (double *)malloc((Nr+1)*(Ne+1)*sizeof(double));
  qnur = (double *)malloc((Nr+1)*(Ne+1)*sizeof(double));
  R_pass = (double *)malloc((Nr+1)*(Ne+1)*sizeof(double));
  y_pass = (double *)malloc((Nth_obs+1)*(Ne+1)*sizeof(double));
  l_pass = (double *)malloc((Nth_obs+1)*(Ne+1)*sizeof(double));
  I_r = (double *)malloc((Nr+1)*sizeof(double));
  I_rph = (double *)malloc((Nr+1)*(Nph+1)*sizeof(double));
  B_factr = (double *)malloc((Nr+1)*sizeof(double));
  G_factr = (double *)malloc((Nr+1)*sizeof(double));
  L_factr = (double *)malloc((Nr+1)*sizeof(double));
  T_factr = (double *)malloc((Nr+1)*sizeof(double));
  I_pass = (double *)malloc((Nr+1)*sizeof(double));
  x_clumps = (double *)malloc((Nclumps+1)*3*sizeof(double));
  r_clumps = (double *)malloc((Nclumps+1)*sizeof(double));
  rho_ijk = (double *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(double));
  T_ijk = (double *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(double));
  bb_ijk = (double *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(double));
  tau_ijk = (double *)malloc((Nr+1)*(Nth+2)*(Nph+1)*sizeof(double));
  ut_ijk = (double *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(double));
  ur_ijk = (double *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(double));
  uz_ijk = (double *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(double));
  up_ijk = (double *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(double));
  corpow_ijk = (double *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(double));
  diskbody_ik = (int *)malloc((Nr+1)*(Nph+1)*sizeof(int));
  sigtau_ik = (double *)malloc((Nr+1)*(Nph+1)*sizeof(double));
  Tdisk_ik = (double *)malloc((Nr+1)*(Nph+1)*sizeof(double));
  emtop_ik = (double *)malloc((Nr+1)*(Nph+1)*sizeof(double));
  embot_ik = (double *)malloc((Nr+1)*(Nph+1)*sizeof(double));
  reftop_ik = (double *)malloc((Nr+1)*(Nph+1)*sizeof(double));
  refbot_ik = (double *)malloc((Nr+1)*(Nph+1)*sizeof(double));
  emtop_elf_ik = (double *)malloc((Nr+1)*(Nph+1)*4*4*sizeof(double));
  embot_elf_ik = (double *)malloc((Nr+1)*(Nph+1)*4*4*sizeof(double));
  reftop_elf_ik = (double *)malloc((Nr+1)*(Nph+1)*4*4*sizeof(double));
  refbot_elf_ik = (double *)malloc((Nr+1)*(Nph+1)*4*4*sizeof(double));
  Gtop_ik = (double *)malloc((Nr+1)*(Nph+1)*sizeof(double));
  Gbot_ik = (double *)malloc((Nr+1)*(Nph+1)*sizeof(double));
  if (imageprint == 1) {
    image = (double *)malloc((Nth_obs+1)*(Ni+1)*(Ni+1)*sizeof(double));
    imagex = (double *)malloc((Nth_obs+1)*(Ni+1)*(Ni+1)*sizeof(double));
    imagey = (double *)malloc((Nth_obs+1)*(Ni+1)*(Ni+1)*sizeof(double));
    spcimage = 
      (double *)malloc((Nth_obs+1)*(Ni+1)*(Ni+1)*(Ne_obs+1)*sizeof(double));
    spcimagex = 
      (double *)malloc((Nth_obs+1)*(Ni+1)*(Ni+1)*(Ne_obs+1)*sizeof(double));
    spcimagey = 
      (double *)malloc((Nth_obs+1)*(Ni+1)*(Ni+1)*(Ne_obs+1)*sizeof(double));
    phimage = 
      (double *)malloc((Nth_obs+1)*(Nph_obs+1)*(Ni+1)*(Ni+1)*sizeof(double));
    phimagex = 
      (double *)malloc((Nth_obs+1)*(Nph_obs+1)*(Ni+1)*(Ni+1)*sizeof(double));
    phimagey = 
      (double *)malloc((Nth_obs+1)*(Nph_obs+1)*(Ni+1)*(Ni+1)*sizeof(double));
    for (it=0;it<=Nth_obs;it++) {
      for (ix=0;ix<=Ni;ix++) {
	for (iy=0;iy<=Ni;iy++) {
	  image[indexi(it,ix,iy)]=0;
	  imagex[indexi(it,ix,iy)]=0;
	  imagey[indexi(it,ix,iy)]=0;
	  for (je=0;je<=Ne_obs;je++) {
	    spcimage[indexspci(it,ix,iy,je)]=0;
	    spcimagex[indexspci(it,ix,iy,je)]=0;
	    spcimagey[indexspci(it,ix,iy,je)]=0;
	  }
	  for (jph=0;jph<=Nph_obs;jph++) {
	    phimage[indexphi(it,jph,ix,iy)]=0;
	    phimagex[indexphi(it,jph,ix,iy)]=0;
	    phimagey[indexphi(it,jph,ix,iy)]=0;
	  }
	}
      }
    }
  }
  for (it=0;it<=Nth_obs;it++) {
    for (je=0;je<=Ne;je++) {
      Ispec[index2(it,je)]=0;
      Lspec[index2(it,je)]=0;
      Qspec[index2(it,je)]=0;
      Uspec[index2(it,je)]=0;
      for (ir=0;ir<=Nr;ir++) {
	Ispecr[indexrth(ir,it,je)]=0;
	Rspecr[indexrth(ir,it,je)]=0;
	Qspecr[indexrth(ir,it,je)]=0;
	Uspecr[indexrth(ir,it,je)]=0;
      }
    }
    for (jph=0;jph<=Nph_obs;jph++) {
      for (jt=0;jt<=Nt;jt++) {
	Iobs[indexph(it,jph,jt)]=0;
      }
      for (je=0;je<=Ne;je++) {
	Ispecp[indexph(it,jph,je)]=0;
      }
    }
  }
  for (ir=0;ir<=Nr;ir++) {
    for (je=0;je<=Ne;je++) {
      Rspec[indexre(ir,je)]=0;
      Cspecr[indexre(ir,je)]=0;
      Inur[indexre(ir,je)]=0;
      qnur[indexre(ir,je)]=0;
    }
  }
  for (ir=0;ir<=Nr;ir++) {
    I_r[ir]=0;
    for (j=0;j<=Nph;j++) {
      I_rph[indexr(ir,j)]=0;
    }
  }
  for (ir=0;ir<=Nr;ir++) {
    for (ith=0;ith<=Nth;ith++) {
      for (iph=0;iph<=Nph;iph++) {
	corpow_ijk[indexijk(ir,ith,iph)]=0;
      }
    }
  }

  a2 = aa*aa;
  Rhor = M+sqrt(M*M-a2);
  kapp_es = 0.4;
  eps_th = 1e-5;
  th0 = 0.08061;
  f_Fe = 0.5; //Fe Kalpha flourescent yield
  T_e0 = 0.25;
  T_e0 = T_cor/1000.; //corona temp in MeV
  T_d0 = 1.0;
  fcseed = 0.1; //fraction of coronal seeds launched
  f_hard = 1.8;
  if (OPT_THIN == 1) f_hard = 1.0;
  if (em_model == 1) spec_model = 1;  //1 for linear energy scale, 2 for log scale
  if (em_model > 1) spec_model = 2;  //1 for linear energy scale, 2 for log scale
  rr_model = 2;  //1 for linear spacing of rr, 2 for log scale
  R_g = 1.45e6*(Mstar/3.33);
  pro = 1.0;
  Z1 = 1.+pow((1-(aa*aa)/(M*M)),1./3.)*
    (pow((1+aa/M),1./3.)+pow((1-aa/M),1./3.));
  Z2 = sqrt(3.*(aa*aa)/(M*M)+Z1*Z1);
  Risco = M*(3.+Z2-pro*sqrt((3.-Z1)*(3.+Z1+2.*Z2)));
  //  printf("a/M = %8.4e\n",aa);
  //  printf("Risco = %12.6e\n",Risco);
  Eisco = (Risco*Risco-2*M*Risco+pro*aa*sqrt(M*Risco))/
    (Risco*sqrt(Risco*Risco-3*M*Risco+pro*2*aa*sqrt(M*Risco)));
  Lisco = pro*sqrt(M*Risco)*(Risco*Risco-pro*2*aa*sqrt(M*Risco)+a2)/
    (Risco*sqrt(Risco*Risco-3*M*Risco+pro*2*aa*sqrt(M*Risco)));
  r_es = Risco*R_g;

  //linear scale
  if (spec_model == 1) {
    e_min = 0.0;
    e_max = 2.0;
    nu0[0]=e_min;
    for (j=1;j<=Ne;j++) {
      nu0[j]=e_min+((double)j)/Ne*(e_max-e_min);
    }
    for (j=1;j<Ne;j++) dnu0[j]=0.5*(nu0[j+1]-nu0[j-1]);
    dnu0[0]=nu0[1]-nu0[0];
    dnu0[Ne]=nu0[Ne]-nu0[Ne-1];
    nui0[0]=e_min;
    for (j=1;j<=Ne_obs;j++) {
      nui0[j]=e_min+((double)j)/Ne_obs*(e_max-e_min);
    }
    for (j=1;j<Ne_obs;j++) dnui0[j]=0.5*(nui0[j+1]-nui0[j-1]);
    dnui0[0]=nui0[1]-nui0[0];
    dnui0[Ne_obs]=nui0[Ne_obs]-nui0[Ne_obs-1];
  }
  //log scale
  if (spec_model == 2) {
    e_min = 0.001;
    if (Mstar > 100) e_min = 0.001; //AGN temp scale
    e_max = 1000.0;
    //e_min = 1e-7;                   //synchrotron seeds
    //e_min = 1e-3;
    e_max = 10000.0;                //hot corona scale
    nu0[0]=e_min;
    for (j=1;j<=Ne;j++) {
      nu0[j]=e_min*pow(10.,((double)j)/Ne*log10(e_max/e_min));
      T_tab[j]=nu0[j];
    }
    for (j=1;j<Ne;j++) dnu0[j]=0.5*(nu0[j+1]-nu0[j-1]);
    dnu0[0]=nu0[1]-nu0[0];
    dnu0[Ne]=nu0[Ne]-nu0[Ne-1];

    nui0[0]=e_min;
    for (j=1;j<=Ne_obs;j++) {
      nui0[j]=e_min*pow(10.,((double)j)/Ne_obs*log10(e_max/e_min));
    }
    for (j=1;j<Ne_obs;j++) dnui0[j]=0.5*(nui0[j+1]-nui0[j-1]);
    dnui0[0]=nui0[1]-nui0[0];
    dnui0[Ne_obs]=nui0[Ne_obs]-nui0[Ne_obs-1];
  }
  calc_Pnorm_T(Pnorm,T_tab,Ne+1);
  Carter = 0;
  t = 0.0;
  dt = 1.0;
  V_tot = 0;
  srand(time(NULL));
  //srand(RUN_ID);
  fov = 20;
  Rshell = 10000;
  Rin = Risco;
  //Rin = 5.; /* Schematic plot */
  Rin = Rhor*1.01;
  Rout_flux = 1000.0;
  //Rout = 200.0;
  //Rin = Risco;
  //Rout = Rin+0.01;
  //Rout = 5.1;
  Redge = Risco;
  if ((em_model == 1.5)||(em_model == 3.5)) Redge = Rin;
  Tmin = Rshell-2.*Rout;
  Tmax = Rshell+100.;
  dTbin = (Tmax-Tmin)/(Nt+1.);
  irstart = 0;
  view_jph = 30;
  //dph = 2.*PI/(N);
  Nbins = (N+1)*(N+1)+1;
  E0 = 3.0;
  z_hat[0]=0;
  z_hat[1]=0;
  z_hat[2]=1;
  get_harm3d_data(rr,tt,pp,rho_ijk,T_ijk,bb_ijk,tau_ijk,ut_ijk,ur_ijk,uz_ijk,up_ijk,
  	  diskbody_ik,sigtau_ik,Tdisk_ik,emtop_ik,embot_ik,reftop_ik,refbot_ik);
  printf("check\n");
  /*
  T_e0 = 1e7;
  T_e = 1e7; 
  for (ir=0;ir<=Nr;ir++) {
    for (ith=0;ith<=Nth;ith++) {
      for (iph=0;iph<=Nph;iph++) {
	if (T_ijk[indexijk(ir,ith,iph)] > T_e) {
	  T_e = T_ijk[indexijk(ir,ith,iph)];
	}
	if ((T_ijk[indexijk(ir,ith,iph)] < T_e0)&&
	    (T_ijk[indexijk(ir,ith,iph)] > 0)) {
	  T_e0 = T_ijk[indexijk(ir,ith,iph)];
	}
      }
    }
  }
  printf("temps %g %g\n", T_e0, T_e);
  */
  nt_spectrum(Risco,rr,drr,nu0,Inur_NT,ntT,qnur);

  //set the radial values of the boundaries of the cells
  rrb[0]=rr[0]-0.5*(rr[1]-rr[0]);
  rrb[Nr+1]=rr[Nr]+0.5*(rr[Nr]-rr[Nr-1]);
  Rin = rrb[0];
  Rout = rrb[Nr+1];
  gamma = pow((Rout/Rin),1./(Nr+1.));
  irstart = 0;
  for (ir=0;ir<=Nr;ir++) {
    rrb[ir+1]=rrb[ir]*gamma;
    drr[ir]=rrb[ir+1]-rrb[ir];
    if (rr[ir] < Rhor) irstart++;
  }
  dph = (pp[Nph]-pp[0])/Nph;
  dth = (tt[Nth]-tt[0])/Nth;
  if (irstart == 0) irstart = 1;
  //outfile = fopen("ph_traj.dat","w");
  surface_tetrads(rr,tt,pp,ut_ijk,ur_ijk,uz_ijk,up_ijk,
		  diskbody_ik,emtop_ik,embot_ik,reftop_ik,refbot_ik,
		  emtop_elf_ik,embot_elf_ik,reftop_elf_ik,refbot_elf_ik,
		  rho_ijk,T_ijk,bb_ijk,Gtop_ik,Gbot_ik,irstart,eps_th);
  
  //irstart = 10;
  //irstop = Nr;
  irstop = irstart;
  for (ir=irstop;ir<=Nr;ir++) {
    if (rr[ir] < Rout_flux) irstop++;
  }
  irstop=irstop-1;
  printf("%ld %g\n", irstop, rr[irstop]);
  //irstart = 100;
  //irstop = 100;
  irstep = 1;
  iphstart = 0;
  iphstop = Nph;
  iphstep = 1;

  flux2 = 0;
  //for (ibottom=0;ibottom<=1;ibottom++) {
  for (ibottom=0;ibottom<=TWO_SIDED;ibottom++) {
  for (ir=irstart;ir<=irstop;ir+=irstep) {
    for (iph=iphstart;iph<=iphstop;iph+=iphstep) {
      srand(ibottom*(Nr+1)*(Nph+1)+ir*(Nph+1)+iph+RUN_ID*0);
      y0[0]=0;
      y0[1]=rr[ir];
      if (ibottom == 0) y0[2]=emtop_ik[indexr(ir,iph)]-eps_th;
      if (ibottom == 1) y0[2]=embot_ik[indexr(ir,iph)]+eps_th;
      y0[3]=pp[iph];
      //lambda = (double)rand()/(RAND_MAX);
      //y0[3]=y0[3]+((int)(lambda*4.))*PI/2.;
      r = y0[1];
      t = y0[2];
      f = y0[3];
      r2 = r*r;
      r3 = r*r2;
      Sig = r*r+a2*cos(t)*cos(t);      //rho^2 in some texts
      Delta = r*r-2*M*r+a2;
      alpha = sqrt(Sig*Delta/(Sig*Delta+2*M*r*(a2+r2)));
      omg = 2.*M*r*aa/(Sig*Delta+2.*M*r*(a2+r2));
      omgb = sqrt((Sig*Delta+2.*M*r*(a2+r2))/Sig*sin(t)*sin(t));
      
      for (j=0;j<=3;j++) part_x0[j]=y0[j];
      calc_g(g_dn,g_up,part_x0);
      if (ibottom == 0) {
	G_fact = Gtop_ik[indexr(ir,iph)];
	for (i=0;i<=3;i++) {
	  for (j=0;j<=3;j++) {
	    e_lf[i][j]=emtop_elf_ik[indexelf(ir,iph,i,j)];
	  }
	}
      }
      if (ibottom == 1) {
	G_fact = Gbot_ik[indexr(ir,iph)];
	for (i=0;i<=3;i++) {
	  for (j=0;j<=3;j++) {
	    e_lf[i][j]=embot_elf_ik[indexelf(ir,iph,i,j)];
	  }
	}
      }
      for (j=0;j<=3;j++) part_v[j]=e_lf[0][j];
      dx_r[0] = 0;
      dx_r[1] = drr[ir];
      dx_r[2] = 0;
      dx_r[3] = 0;
      dx_th[0] = 0;
      dx_th[1] = 0;
      dx_th[2] = dth;
      dx_th[3] = 0;
      dx_p[0] = 0;
      dx_p[1] = 0;
      dx_p[2] = 0;
      dx_p[3] = dph;

      if (ibottom <= 1) {
	T_e0 = Tdisk_ik[indexr(ir,iph)];
	//T_e0 = ntT[ir];
	//printf("%d %d %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %d\n",
	//     ir,iph,rr[ir],y0[2],y0[3],part_v[0],part_v[1],
	//     part_v[2],part_v[3],T_e0,diskbody_ik[indexr(ir,iph)]);
	
	T_e0=T_e0*f_hard;
	for (j=0;j<=Ne;j++) {
	  if (T_e0 <= 0) {
	    Inur[indexre(ir,j)]=0;
	  }
	  if (T_e0 > 0) {
	    x = 1.e3*nu0[j]/(kB_ev*T_e0);
	    Inur[indexre(ir,j)]=pow(nu0[j],3.)/(exp(x)-1.)/pow(f_hard,4.);
	    //kap_ff = 1.5e25*rhoc_r[i]*pow(Ts_r[i],-3.5)*pow(x,-3.)*(1.-exp(-x));
	    //qnur[indexre(i,j)]=kap_es/(kap_es+kap_ff);
	  }
	}
	flux = 0;
	for (j=0;j<=Ne;j++) flux+=Inur[indexre(ir,j)]*dnu0[j];
	//printf("%12.5e %12.5e\n",Inur[indexre(ir,50)],flux*PI);
	if (fmod(ir,10) < 0) {
	  printf("e_tlf: %12.5e %12.5e %12.5e %12.5e\n",
		 e_lf[0][0],e_lf[0][1],e_lf[0][2],e_lf[0][3]);
	  printf("e_xlf: %12.5e %12.5e %12.5e %12.5e\n",
		 e_lf[1][0],e_lf[1][1],e_lf[1][2],e_lf[1][3]);
	  printf("e_ylf: %12.5e %12.5e %12.5e %12.5e\n",
		 e_lf[2][0],e_lf[2][1],e_lf[2][2],e_lf[2][3]);
	  printf("e_zlf: %12.5e %12.5e %12.5e %12.5e\n",
		 e_lf[3][0],e_lf[3][1],e_lf[3][2],e_lf[3][3]);
	  printf("dx_r : %12.5e\n", 0.5*(emtop_ik[indexr(ir+1,iph)]-
					 emtop_ik[indexr(ir-1,iph)]));
	  printf("dx_p : %12.5e\n", 0.5*(emtop_ik[indexr(ir,iph+1)]-
					 emtop_ik[indexr(ir,iph-1)]));
	  
	}
	G_factr[ir] = G_fact;
	T_fact = 1./part_v[0];
	T_factr[ir] = T_fact;
      }

      //printf("%g %g %g\n",e_x_hat[0],e_x_hat[1],e_x_hat[2]);
      //printf("%g %g %g\n",e_y_hat[0],e_y_hat[1],e_y_hat[2]);
      //printf("%g %g %g %g\n",part_v[0],part_v[1],part_v[2],part_v[3]);
      //printf("%g %g %g %g\n",part_v_hat[0],part_v_hat[1],part_v_hat[2],part_v_hat[3]);

      //printf("%12.5e %12.5e %12.5e %12.5e\n",
      //     r,beta,G_fact,T_fact);
      ith = 0;
      phi = 0.;
      tau_tt = 0.;
      Nescape = 0.;
      Ndisk = 0.;
      Nscat = 0.;
      Ncapt = 0.;

      for (it=0;it<=N;it++) {
	for (ip=0;ip<=N;ip++) {
	  ith++;
	  raydex = ith;
	  moddex = raydex%numprocs;
	  err = erro;
	  dt = 1.0;
	  l_tot=0;
	  Fe_phot = 0;
	  tau_tot = 0;
	  steps = 0;
	  dscat = 0;
	  iscat = 0;
	  just_scattered = 0;

	  //LAUNCH DISK PHOTON
	  if ((ibottom == 0)||(ibottom == 1)) {
	    lambda = (double)rand()/(RAND_MAX);
	    cth = lambda; //only consider photons moving in +z direction
	    //if ((rr[ir]<Redge)&&(atm_model==1.5)) cth = 1.-2.*lambda;
	    sth = sqrt(1.-cth*cth);
	    lambda = (double)rand()/(RAND_MAX);	
	    phi = 2.*PI*lambda;
	    
	    r = y0[1];
	    t = y0[2];
	    f = y0[3];

	    p_hat[0] = sth*cos(phi);
	    p_hat[1] = sth*sin(phi);
	    p_hat[2] = cth;

	    cross(p_hat,z_hat,f_hat);
	    normalize(f_hat);
	    //calc_g(g_dn,g_up,part_x0);
	    cth0 = fabs(p_hat[2]);
    
	    //Cartesian coordinate area
	    //A_fact = dph*drr[ir];
	    //These factors should now be included in G_fact
	    //Factor of 8 goes from PI/2 half-plane to full 2PI disk, both sides
	    if (TWO_SIDED == 0) A_fact = 8.0;
	    if (TWO_SIDED > 0) A_fact = 4.0;
	    //cgs area conversion
	    A_fact = A_fact*R_g*R_g;
	    //cos theta factor for optically thick planar emitter
	    A_fact = A_fact*cth0;
	    //dOmega factor for solid angle of each ray
	    A_fact = A_fact*(2*PI)/((double)(N+1)*(N+1));
	    
	    //CALCULATE INITIAL POLARIZATION AND LIMB-DARKENING FROM CHANDRA (1960)
	    chandra_limit(cth0, &deg, &D_fact);

	    //INITIALY UNPOLARIZED LIGHT FROM DISK
	    //deg = 0;
	    //CORRECTION FOR TAU=1
	    //deg = deg/1.25;
	    //D_fact = D_fact*(0.85+0.3*cth0);
	    A_fact = A_fact*D_fact;
	    
	    //isotropic optically thin
	    if (OPT_THIN == 1) A_fact = A_fact/cth0/2./D_fact;
	    //isotropic optically thick:
	    //A_fact = A_fact/D_fact;
	    
	    //Proper area in emitter frame
	    A_fact = A_fact*G_factr[ir];
	    
	    //Time dilation: observed intensity should be lower by 
	    //a factor of dtau/dt
	    A_fact = A_fact*T_factr[ir];
	    
	    //Bnu_ is spectral brightness measured by emitter
	    for (iv = 0;iv<=Ne;iv++) Bnu_[iv]=0;
	    
	    //LINE EMISSION
	    if (em_model == 1) {
	      E_i = 1.0;
	      for (iv=0;iv<=Ne;iv++) {
		nu[iv]=nu0[iv];
		dnu[iv]=dnu0[iv];
		Bnu_[iv]=0;
	      }
	      //uniform emission
	      //B_fact = 1.0;
	      //1/r^2 emission
	      B_fact = 1.0/(rr[ir]*rr[ir]);
	      iv=(int)floor((E_i-e_min)/(e_max-e_min)*Ne);
	      Bnu_[iv]=B_fact;
	    }
	    
	    //THERMAL EMISSION
	    if (em_model >= 2) {
	      for (iv=0;iv<=Ne;iv++) {
		nu[iv] = nu0[iv];
		dnu[iv] = dnu0[iv];
		//Brightness in units of [ergs/cm^2/s/Hz/Sr]	      
		Bnu_[iv]=2.*1.04e5*Inur[indexre(ir,iv)];
	      }
	      if ((it < 0)&&(ip < 0)) {
		printf("%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
		       nu[40],nu[50],nu[60],nu[70],nu[80],nu[90],nu[100]);
		printf("%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
		       Bnu_[40],Bnu_[50],Bnu_[60],Bnu_[70],Bnu_[80],Bnu_[90],Bnu_[100],
		       T_e0);
	      }
	      if (diskbody_ik[indexr(ir,iph)] == 0) A_fact =0;
	      if (Tdisk_ik[indexr(ir,iph)] == 0) A_fact = 0;
	    }
	  }
	  //LAUNCH CORONA PHOTON
	  if ((ibottom == 2)||(ibottom == 3)) {
	    lambda = (double)rand()/(RAND_MAX);
	    cth = 1.-2.*lambda; //photons emitted isotropically
	    sth = sqrt(1.-cth*cth);
	    lambda = (double)rand()/(RAND_MAX);	
	    phi = 2.*PI*lambda;
	    
	    lambda = (double)rand()/(RAND_MAX);	
	    //upper corona
	    if (ibottom == 2) { 
	      dcth = (emtop_ik[indexr(ir,iph)]-th0)/((double)(N+1)*(N+1));
	      y0[2] = th0+((double)it*(N+1)+(double)ip+1.)*dcth;
	    }
	    //lower corona
	    if (ibottom == 3) { 
	      dcth = (PI-th0-embot_ik[indexr(ir,iph)])/((double)(N+1)*(N+1));
	      y0[2] = PI-th0-((double)it*(N+1)+(double)ip+1.)*dcth;
	    }
	    r = y0[1];
	    t = y0[2];
	    f = y0[3];
	    for (j=0;j<=3;j++) part_x[j]=y0[j];
	    calc_g(g_dn_ph,g_up_ph,part_x);
	    lookup_data(part_x,rr,tt,pp,g_dn_ph,
			rho_ijk,T_ijk,bb_ijk,ut_ijk,ur_ijk,uz_ijk,up_ijk,
			weights,&rho,&T_e,&Bmag,part_v);
	    //printf("data %g %g %g\n",y0[2],T_e,part_v[0]);
	    part_p[0] = g_dn_ph[0][0]*part_v[0]+g_dn_ph[0][3]*part_v[3];
	    part_p[1] = g_dn_ph[1][1]*part_v[1];
	    part_p[2] = g_dn_ph[2][2]*part_v[2];
	    part_p[3] = g_dn_ph[3][0]*part_v[0]+g_dn_ph[3][3]*part_v[3];
	    T_fact = 1./part_v[0];

	    p_hat[0] = sth*cos(phi);
	    p_hat[1] = sth*sin(phi);
	    p_hat[2] = cth;

	    cross(p_hat,z_hat,f_hat);
	    normalize(f_hat);
	    cth0 = fabs(p_hat[2]);
    
	    //Factor of 4 goes from PI/2 half-plane to full 2PI disk
	    A_fact = 4.0;
	    //cgs area conversion
	    A_fact = A_fact*R_g*R_g;
	    //dOmega factor for solid angle of each ray
	    A_fact = A_fact*4*PI;
	    //factor from sampling entire corona volume in (N+1)^2 slices
	    A_fact = A_fact*(dcth/dth);
	    //factor for sparce sampling of coronal photons
	    A_fact = A_fact/fcseed;
	    lambda = (double)rand()/(RAND_MAX);
	    if (lambda < fcseed) A_fact = 0;
	    deg = 0;
	    
	    //Proper volume in emitter frame
	    calc_dV(e_lf,w_lf,g_up_ph,g_dn_ph,part_p,part_v,
		    dx_r,dx_th,dx_p,&V_fact,&G_fact);
	    //For consistency with disk surface, treat coronal volume
	    //element as a surface with area G_fact and thickness h_fact,
	    //such that V_fact = G_fact*h_fact
	    h_fact = V_fact/G_fact*1.475e5*3.*Mstar;

	    dl = pow(V_fact,1./3.);

	    //printf("%g %g %g %g %g %g\n", r, embot_ik[indexr(ir,iph)],t, V_fact, dth, dcth);
	    A_fact = A_fact*dl*dl;
	    V_tot += V_fact*4.*(dcth/dth);

	    //Time dilation: observed intensity should be lower by 
	    //a factor of dtau/dt
	    A_fact = A_fact*T_fact;
	    
	    //Bnu_ is spectral brightness measured by emitter
	    for (iv = 0;iv<=Ne;iv++) Bnu_[iv]=0;
	    
	    jth = (int)(fabs(cos(y0[2]))*(Nth+1.));
	    if (TWO_SIDED >= 1) {
	      jth = (int)((cos(y0[2])+1.)/2.*(Nth+1.));
	    }
	    if (jth > Nth) jth = Nth;
	    
	    //THERMAL CYCLOTRON EMISSION
	    if (em_model >= 20) {
	      gam2 = pow((kB_ev*T_e/5.11e5+1.),2.);
	      nuc = qe*sqrt(Bmag)/(2*PI*me*cc); //nu_c in Hz
	      nuc = nuc/2.42e17; //nu_c in keV
	      if (nuc == 0) nuc = 1.e-5;
	      n_e = rho/mp;
	      P0 = lookup_Pnorm(Pnorm,nu0,T_e*kB_ev*1e-3,Ne);
	      printf("%g %g %g\n",Bmag,gam2,nuc);
	      if (P0 > 0) P0 = 4./3.*sigma_T*cc*(gam2-1)*(Bmag/(8.*PI))*n_e/
		(nuc*2.42e17*P0*4.*PI);
	      Iptot = 0;
	      for (iv=0;iv<=Ne;iv++) {
		nu[iv] = nu0[iv];
		dnu[iv] = dnu0[iv];
		//blackbody brightness in units of [ergs/cm^2/s/Hz/Sr]	      
		B_nu[iv] = 0.;
		if (T_e > 0) {
		  x = nu[iv]/(kB_ev*T_e*1e-3);
		  B_nu[iv]=2.*1.04e5*(nu[iv]*nu[iv]*nu[iv])/(exp(x)-1.);
		}
		//cyclotron emissivity in units of [ergs/cm^3/s/Hz/Sr]
		j_nu[iv]= 0.;
		x = nu[iv]/nuc;
		if ((T_e > 1e7)&&(Bmag > 0)&&(x <= 1)) 
		  j_nu[iv]=P0*F_cycl(x,kB_ev*T_e*1e-3);
		//cyclotron self-absorption in units of [1/cm]
		alpha_nu[iv] = 0;
		if (B_nu[iv] > 0) alpha_nu[iv]=j_nu[iv]/B_nu[iv];
		//Brightness escaping the corona cell
		Bnu_[iv]=B_nu[iv]*(1.-exp(-alpha_nu[iv]*dl*R_g)); 
		//optically thin cyclotron
		Bnu_[iv]=j_nu[iv]*dl*R_g; 
		for (i=0;i<=1;i++) {
		  for (j=0;j<=1;j++) {
		    for (k=0;k<=1;k++) {
		      wght = weights[i+6]*weights[j+8]*weights[k+10];
		      corpow_ijk[indexijk((int)weights[i],(int)weights[j+2],
					  (int)weights[k+4])] 
			+= wght*Bnu_[iv]*dnu0[iv]*A_fact;
		    }
		  }
		}
		Iptot+=dnu0[iv]*j_nu[iv]*2.42e17*4.*PI;
	      }
	      if ((it < 0)&&(ip < 0)) {
		printf("%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
		       nu[10],nu[20],nu[30],nu[70],nu[80],nu[90],nu[100]);
		printf("%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
		       Bnu_[10],Bnu_[20],Bnu_[30],Bnu_[70],Bnu_[80],Bnu_[90],Bnu_[100],
		       wc,gam2,T_e,sqrt(Bmag));
	      }
	    }
	    //THERMAL SYNCHROTRON EMISSION
	    if (em_model >= 20) {
	      gam2 = pow((kB_ev*T_e/5.11e5+1.),2.);
	      wc = 3.*PI/8.*gam2*qe*sqrt(Bmag)/(2.*PI*me*cc); //omega_c in Hz
	      wc = wc/2.42e17; //omega_c in keV
	      n_e = rho/mp;
	      //printf("%g %g %g\n",Bmag,gam2,wc);
	      P0 = 4./sqrt(3.)/(PI*PI)*(qe*qe*qe/(me*cc*cc))* 
		sqrt(Bmag)*(1.-1./gam2)*n_e/(4.*PI);
	      if (wc == 0) wc = 1.e-5;
	      Iptot = 0;
	      for (iv=0;iv<=Ne;iv++) {
		nu[iv] = nu0[iv];
		dnu[iv] = dnu0[iv];
		//blackbody brightness in units of [ergs/cm^2/s/Hz/Sr]	      
		B_nu[iv] = 0.;
		if (T_e > 0) {
		  x = nu[iv]/(kB_ev*T_e*1e-3);
		  B_nu[iv]=2.*1.04e5*(nu[iv]*nu[iv]*nu[iv])/(exp(x)-1.);
		}
		//synchrotron emissivity in units of [ergs/cm^3/s/Hz/Sr]
		j_nu[iv]= 0.;
		x = nu[iv]/wc;
		if ((T_e > 0)&&(Bmag > 0)) 
		  j_nu[iv]=P0*F_sync(x);
		//synchrotron self-absorption in units of [1/cm]
		alpha_nu[iv] = 0;
		if (B_nu[iv] > 0) alpha_nu[iv]=j_nu[iv]/B_nu[iv];
		//Brightness escaping the corona cell
		Bnu_[iv]=B_nu[iv]*(1.-exp(-alpha_nu[iv]*dl*R_g)); 
		//optically thin synchrotron
		Bnu_[iv]=j_nu[iv]*dl*R_g; 
		for (i=0;i<=1;i++) {
		  for (j=0;j<=1;j++) {
		    for (k=0;k<=1;k++) {
		      wght = weights[i+6]*weights[j+8]*weights[k+10];
		      corpow_ijk[indexijk((int)weights[i],(int)weights[j+2],
					  (int)weights[k+4])] 
			+= wght*Bnu_[iv]*dnu0[iv]*A_fact;
		    }
		  }
		}

		//corpow_ijk[indexijk(ir,jth,iph)]+=Bnu_[iv]*dnu0[iv]*A_fact;
		Iptot+=dnu0[iv]*P0*n_e*F_sync(x)*2.42e17;
	      }
	      //printf("%10.5e\n",Iptot/((4./3.*(gam2-1.))*cc*6.66e-25*Bmag/(8.*PI)*n_e));
	      if ((it < 0)&&(ip < 0)) {
		printf("%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
		       nu[10],nu[20],nu[30],nu[70],nu[80],nu[90],nu[100]);
		printf("%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
		       Bnu_[10],Bnu_[20],Bnu_[30],Bnu_[70],Bnu_[80],Bnu_[90],Bnu_[100],
		       wc,gam2,T_e,sqrt(Bmag));
	      }
	    }
	    //THERMAL BREMSSTRAHLUNG EMISSION
	    if (em_model >= 20) {
	      for (iv=0;iv<=Ne;iv++) {
		nu[iv] = nu0[iv];
		dnu[iv] = dnu0[iv];
		//Bnu_[iv] = 0;
	      }
	      if (T_e >0) {
		Iptot = 0;
		for (iv=0;iv<=Ne;iv++) {
		  //Spectral emissivity in units of [ergs/cm^3/s/Hz]
		  x = nu[iv]/(kB_ev*T_e*1e-3);
		  n_e = rho/mp;
		  j_nu[iv]=6.8e-38*n_e*n_e/sqrt(T_e)*exp(-x)/(4*PI);
		  Bnu_[iv]+=j_nu[iv]*dl*R_g; 
		  Iptot+=Bnu_[iv]*dnu0[iv];
		  for (i=0;i<=1;i++) {
		    for (j=0;j<=1;j++) {
		      for (k=0;k<=1;k++) {
			wght = weights[i+6]*weights[j+8]*weights[k+10];
			corpow_ijk[indexijk((int)weights[i],(int)weights[j+2],
					    (int)weights[k+4])] 
			  += wght*Bnu_[iv]*dnu0[iv]*A_fact;
		      }
		    }
		  }
		  if ((it < 0)&&(ip < 0)) {
		    printf("%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
			   nu[40],nu[50],nu[60],nu[70],nu[80],nu[90],nu[100]);
		    printf("%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
			   Bnu_[40],Bnu_[50],Bnu_[60],Bnu_[70],Bnu_[80],Bnu_[90],Bnu_[100],
			   n_e,T_e);
		  }
		}
	      }
	    }
	    //if (Iptot > 10) printf("launch %g %g\n", t, Iptot);

	    for (i=0;i<=3;i++) {
	      for (j=0;j<=3;j++) {
		g_dn[i][j]=g_dn_ph[i][j];
		g_up[i][j]=g_up_ph[i][j];
	      }
	    }
	    //printf("%12.5e %12.5e %12.5e %12.5e %12.5e\n",
	    //   A_fact,cth0,D_fact,G_fact,T_fact);
	  }
	  
	  ph_v[0] = (e_lf[0][0]*(1.)+e_lf[1][0]*p_hat[0]+
		     e_lf[2][0]*p_hat[1]+e_lf[3][0]*p_hat[2]);
	  ph_v[1] = (e_lf[0][1]*(1.)+e_lf[1][1]*p_hat[0]+
		     e_lf[2][1]*p_hat[1]+e_lf[3][1]*p_hat[2]);
	  ph_v[2] = (e_lf[0][2]*(1.)+e_lf[1][2]*p_hat[0]+
		     e_lf[2][2]*p_hat[1]+e_lf[3][2]*p_hat[2]);
	  ph_v[3] = (e_lf[0][3]*(1.)+e_lf[1][3]*p_hat[0]+
		     e_lf[2][3]*p_hat[1]+e_lf[3][3]*p_hat[2]);
	  f0_[0] = (e_lf[0][0]*(0.)+e_lf[1][0]*f_hat[0]+
		    e_lf[2][0]*f_hat[1]+e_lf[3][0]*f_hat[2]);
	  f0_[1] = (e_lf[0][1]*(0.)+e_lf[1][1]*f_hat[0]+
		    e_lf[2][1]*f_hat[1]+e_lf[3][1]*f_hat[2]);
	  f0_[2] = (e_lf[0][2]*(0.)+e_lf[1][2]*f_hat[0]+
		    e_lf[2][2]*f_hat[1]+e_lf[3][2]*f_hat[2]);
	  f0_[3] = (e_lf[0][3]*(0.)+e_lf[1][3]*f_hat[0]+
		    e_lf[2][3]*f_hat[1]+e_lf[3][3]*f_hat[2]);
	  //E0 = dot_g4(g_dn_ph,ph_v,ph_v);
	  //if (fabs(E0) > 1e-3) 
	  //printf("%d %d %d %g %g\n",steps,iscat,dscat,dot_g4(g_dn_ph,ph_v,ph_v),yn[1]);
	  kap10 = aa*cos(t)*((ph_v[0]*f0_[1]-ph_v[1]*f0_[0])
			     +aa*sin(t)*sin(t)*(ph_v[1]*f0_[3]-ph_v[3]*f0_[1]))
	    +r*sin(t)*((r*r+aa*aa)*(ph_v[3]*f0_[2]-ph_v[2]*f0_[3])
		       -aa*(ph_v[0]*f0_[2]-ph_v[2]*f0_[0]));
	  kap20 = r*(ph_v[0]*f0_[1]-ph_v[1]*f0_[0]
		     +aa*sin(t)*sin(t)*(ph_v[1]*f0_[3]-ph_v[3]*f0_[1]))-
	    aa*sin(t)*cos(t)*((r*r+aa*aa)*(ph_v[3]*f0_[2]-ph_v[2]*f0_[3])
			      -aa*(ph_v[0]*f0_[2]-ph_v[2]*f0_[0]));
	  kap1i = kap10;
	  kap2i = kap20;

	  ph_p[0] = g_dn[0][0]*ph_v[0]+g_dn[0][3]*ph_v[3];
	  ph_p[1] = g_dn[1][1]*ph_v[1];
	  ph_p[2] = g_dn[2][2]*ph_v[2];
	  ph_p[3] = g_dn[3][0]*ph_v[0]+g_dn[3][3]*ph_v[3];
	  pdv_em = (ph_p[0]*part_v[0]+ph_p[1]*part_v[1]+
		    ph_p[2]*part_v[2]+ph_p[3]*part_v[3]);

	  for (j=0;j<=3;j++) y0[j+4]=ph_p[j];
	  Carter0 = y0[6]*y0[6]+cos(y0[2])*cos(y0[2])*
	    (y0[7]*y0[7]/(sin(y0[2])*sin(y0[2]))-y0[4]*y0[4]*a2);
	  for (j=0;j<=7;j++) {
	    y[j]=y0[j];
	    y1[j]=y[j];
	    yn[j]=y[j];
	  }
	  calc_g(g_dn_ph,g_up_ph,yn);
	  ph_v[0] = g_up_ph[0][0]*yn[4]+g_up_ph[0][3]*yn[7];
	  ph_v[1] = g_up_ph[1][1]*yn[5];
	  ph_v[2] = g_up_ph[2][2]*yn[6];
	  ph_v[3] = g_up_ph[0][3]*yn[4]+g_up_ph[3][3]*yn[7];
	  for (j=0;j<=3;j++) {
	    p_[j]=yn[j+4];
	  }
	  //printf("%d %d %g %g\n",it,ip,dot_g4(g_dn_ph,part_v,part_v),dot_g4(g_up_ph,p_,p_));

	  //Stop the photon if it crosses the photosphere, the horizon,
	  //escapes the system, or takes too many steps.
	  while ((yn[1]>1.02*Rhor)&&(yn[1]<Rshell*0.99999)&&(steps<1000)&&(A_fact>0)) {
	    for (j=0;j<=7;j++) {
	      y[j]=yn[j];
	    }
	    if (-yn[4]>e_max) A_fact = 0;
	    acc = 0;
	    while ((acc < 1)&&(steps < 1010)) {
	      //calc_e4(e_lf,w_lf,g_up,g_dn,part_p,part_v,dx_r,dx_p,ibottom,&G_fact);
	      if (erro > 0) {
		dt = 0.9*dt*pow((erro/err),0.2);
	      } else {
		dt = 0.9*dt*pow((erro/err),0.25);
	      }
	      if ((yn[1] < rr[Nr/2])&&(dt > yn[1]/20.)) dt = yn[1]/20.;
	      cashkarp(dt,y,y_ck,del_ck);
	      for (j=0;j<=7;j++) {
		del[j]=del_ck[j];
		y2[j]=y_ck[j];
	      }
	      
	      delmag = 0;
	      for (j=0;j<=7;j++) {
		delmag+=del[j]*del[j]/(y2[j]*y2[j]+erro);
	      }
	      err = sqrt(delmag);
	      if (err < erro) {
		acc = 1;
	      }
	      //printf("%d %d %d %g\n",steps,iscat,dscat,err);
	      //if (steps > 600) printf("%d %d %.6g %.6g %6g\n",it,ip,y2[6],dt,
	      //		      (Carter-Carter0)/Carter0);
	      //else printf("%d %d %.6g %.6g\n",it,steps,y2[0],err);
	      if (err < emin) {
		//	    printf("%d %d %.6g %.6g\n",it,ip,y2[0],err);
		err = emin;
	      }
	      steps=steps+1;
	      //if (steps < 2000) fprintf(outfile,"%d %d %12.4e %12.4e %12.4e %12.4e\n",
	      //	       ip,steps,y[1],y[2],y[3],A_fact);
	    }
	    /****************INTERSECTION WITH THE DISK****************/
	    if ((yn[1] < Rout)&&(y2[1] < Rout)&&(y2[1] > Rhor)&&(OPT_THIN == 0)) {
	      just_scattered=0; 
	      //ir_old=(int)(log(yn[1]/rr[0])/log(rr[Nr]/rr[0])*Nr);
	      ir_old=(int)(log(yn[1]/Rin)/log(Rout/Rin)*(Nr+1));
	      if (ir_old > Nr) ir_old = Nr;
	      f_old = yn[3];
	      while (f_old < 0) f_old+=2.*PI;
	      f_old = fmod(f_old,(PI/2.));
	      jph_old =(int)(f_old/dph);
	      if (jph_old > Nph) jph_old = 0;
	      old_zone = 0;
	      if ((cos(yn[2]) < cos(emtop_ik[indexr(ir_old,jph_old)]))&&
		  (cos(yn[2]) > cos(embot_ik[indexr(ir_old,jph_old)])))
		old_zone = 1;
	      if ((cos(yn[2]) < cos(reftop_ik[indexr(ir_old,jph_old)]))&&
		  (cos(yn[2]) > cos(refbot_ik[indexr(ir_old,jph_old)])))
		old_zone = 2;
	      //ir_new=(int)(log(y2[1]/rr[0])/log(rr[Nr]/rr[0])*Nr);
	      ir_new=(int)(log(y2[1]/Rin)/log(Rout/Rin)*(Nr+1));
	      if (ir_new > Nr) ir_new = Nr;
	      f_new = y2[3];
	      while (f_new < 0) f_new+=2.*PI;
	      f_new = fmod(f_new,(PI/2.));
	      jph_new =(int)(f_new/dph);
	      if (jph_new > Nph) jph_new = 0;
	      new_zone = 0;
	      if ((cos(y2[2]) < cos(emtop_ik[indexr(ir_new,jph_new)]))&&
		  (cos(y2[2]) > cos(embot_ik[indexr(ir_new,jph_new)])))
		new_zone = 1;
	      if ((cos(y2[2]) < cos(reftop_ik[indexr(ir_new,jph_new)]))&&
		  (cos(y2[2]) > cos(refbot_ik[indexr(ir_new,jph_new)])))
		new_zone = 2;
	      if ((old_zone > 10)&&(steps <500)) {
		//if ((it == 2)&&(ip == 6)) {
		printf("embedded %ld %d %d %10.4e %10.4e %ld %d %10.4e %10.4e\n",
		       ir,ir_old, jph_old,yn[2],yn[3],steps,jph_new,y2[2],
		       emtop_ik[indexr(ir_old,jph_old)]);
	      }
	      if ((old_zone == 0)&&(new_zone >= 1)) {
		//printf("crossover %d %d %d %d\n",ir_old,jph_old,ir_new,jph_new);
		if (cos(yn[2]) > cos(emtop_ik[indexr(ir_old,jph_old)])) {
		  isbottom = 0;
		}
		if (cos(yn[2]) < cos(embot_ik[indexr(ir_old,jph_old)])) {
		  isbottom = 1;
		}
		for (j=0;j<=7;j++) yn[j]=y2[j];
		r = yn[1];
		r2 = r*r;
		cth = cos(yn[2]);
		sth = sin(yn[2]);
		calc_g(g_dn_ph,g_up_ph,yn);
		ph_v[0] = g_up_ph[0][0]*yn[4]+g_up_ph[0][3]*yn[7];
		ph_v[1] = g_up_ph[1][1]*yn[5];
		ph_v[2] = g_up_ph[2][2]*yn[6];
		ph_v[3] = g_up_ph[0][3]*yn[4]+g_up_ph[3][3]*yn[7];
		for (j=0;j<=3;j++) {
		  part_x[j]=yn[j];
		  k_[j]=ph_v[j];
		  p_[j]=yn[j+4];
		}
		lookup_data(part_x,rr,tt,pp,g_dn_ph,
			    rho_ijk,T_ijk,bb_ijk,ut_ijk,ur_ijk,uz_ijk,up_ijk,
			    weights,&rho,&T_e,&Bmag,v_);
		E0 = g_dn_ph[0][0]*v_[0]*v_[0]+2.*g_dn_ph[0][3]*v_[0]*v_[3]
		  +g_dn_ph[1][1]*v_[1]*v_[1]+g_dn_ph[2][2]*v_[2]*v_[2]
		  +g_dn_ph[3][3]*v_[3]*v_[3];
		if (E0 > -0.99) printf("vdisk %ld %g %g %g\n",steps,E0,r,yn[2]);

		
		po_[0] = g_dn_ph[0][0]*v_[0]+g_dn_ph[0][3]*v_[3];
		po_[1] = g_dn_ph[1][1]*v_[1];
		po_[2] = g_dn_ph[2][2]*v_[2];
		po_[3] = g_dn_ph[0][3]*v_[0]+g_dn_ph[3][3]*v_[3];
		
		//rdsh is redshift of incident photons as measured by disk
		rdsh = (v_[0]*p_[0]+v_[1]*p_[1]+v_[2]*p_[2]+v_[3]*p_[3])/
		  pdv_em;
		for (iv=0;iv<=Ne;iv++) {
		  nu[iv] = nu[iv]*rdsh;
		  dnu[iv] = dnu[iv]*rdsh;
		}
		
		//calculate the incident flux and spectrum on the disk

		//A_fact2 = dph*drr[ir_new];
		//These factors should now be included in G_fact

		//Factor of 8 goes from PI/2 half-plane to full 2PI disk, both sides
		if (TWO_SIDED == 0) A_fact2 = 8.0;
		if (TWO_SIDED > 0) A_fact2 = 4.0;
		//cgs area conversion
		A_fact2 = A_fact2*R_g*R_g;

		if (isbottom == 0) {
		  G_fact = Gtop_ik[indexr(ir_new,jph_new)];
		  T_fact = 1./emtop_elf_ik[indexelf(ir_new,jph_new,0,0)];
		}
		if (isbottom == 1) {
		  G_fact = Gbot_ik[indexr(ir_new,jph_new)];
		  T_fact = 1./embot_elf_ik[indexelf(ir_new,jph_new,0,0)];
		}
		//Proper area in absorber frame
		A_fact2 = A_fact2*G_fact;
		//Time dilation: observed intensity should be lower by 
		//a factor of dtau/dt
		A_fact2 = A_fact2*T_fact;
		
		for (iv=0;iv<Ne;iv++) {
		  I_r[ir_new] += Bnu_[iv]*dnu[iv]*A_fact/A_fact2;
		  I_rph[indexr(ir_new,jph_new)] += 
		    Bnu_[iv]*dnu[iv]*A_fact/A_fact2;
		  //linear energy scale
		  if (spec_model == 1) 
		    je=(int)floor((nu[iv]-e_min)/(e_max-e_min)*Ne);
		  //log energy scale
		  if (spec_model == 2) 
		    je=(int)floor((Ne)*log10(nu[iv]/e_min)/log10(e_max/e_min));
		  if ((je >= 0)&(je <= Ne)) {
		    Rspec[indexre(ir_new,je)] += 
		      Bnu_[iv]*dnu[iv]*A_fact/A_fact2/dnu0[je];	
		  }
		  //if ((it == 10)&&(ip == 10)) 
		  //printf("%d %d %g %g %g %g\n",ir_new,je,yn[1],nu0[iv],Bnu_[iv],
		  //   Rspec[indexre(ir_new,je)]);
		}
		
		//TRANSFORM PHOTON 4-MOMENTUM TO LOCAL TETRAD AND SOLVE FOR F_
		
		bdata[0]=kap1i;
		bdata[1]=kap2i;
		bdata[2]=0;
		bdata[3]=0;
		for (i=0;i<=3;i++) {
		  for (j=0;j<=3;j++) {
		    adata[i*4+j]=0;
		  }
		}
		adata[0*4+0]=-aa*cth*k_[1]+r*sth*aa*k_[2];
		adata[0*4+1]=aa*cth*k_[0]-aa*sth*sth*k_[3];
		adata[0*4+2]=r*sth*((r2+a2)*k_[3]-aa*k_[0]);
		adata[0*4+3]=a2*cth*sth*sth*k_[1]+r*sth*(-(r2+a2)*k_[2]);
		adata[1*4+0]=r*(-k_[1])+a2*sth*cth*k_[2];
		adata[1*4+1]=r*(k_[0]-aa*sth*sth*k_[3]);
		adata[1*4+2]=aa*sth*cth*((r2+a2)*k_[3]-aa*k_[0]);
		adata[1*4+3]=r*aa*sth*sth*k_[1]-aa*sth*cth*(-(r2+a2)*k_[2]);
		for (j=0;j<=3;j++) {
		  adata[2*4+j]=p_[j];
		  adata[3*4+j]=po_[j];
		}
		for (i=0;i<=3;i++) {
		  for (j=0;j<=3;j++) {
		    adata1[i+1][j+1]=adata[i*4+j];
		  }
		  bdata1[i+1]=bdata[i];
		}
		ludcmp_js(adata1,4,indx,&dd);
		lubksb_js(adata1,4,indx,bdata1);
		for (j=0;j<=3;j++) f_[j]=bdata1[j+1];
		
		calc_e2(e_lfs,w_lfs,g_up_ph,g_dn_ph,po_,v_);
		for (i=0;i<=3;i++) {
		  for (j=0;j<=3;j++) {
		    adata1[i+1][j+1]=e_lfs[j][i];
		  }
		}
		ludcmp_js(adata1,4,indx,&dd);
		for (j=0;j<=3;j++) bdata1[j+1]=k_[j];
		lubksb_js(adata1,4,indx,bdata1);
		for (j=0;j<=3;j++) ph_v_hat[j]=bdata1[j+1];
		for (j=0;j<=3;j++) bdata1[j+1]=f_[j];
		lubksb_js(adata1,4,indx,bdata1);
		for (j=0;j<=3;j++) f_v_hat[j]=bdata1[j+1];
		
		E_i = ph_v_hat[0];
		E_i = -yn[4];
		for (j=1;j<=3;j++) {
		  f_v_hat[j]=f_v_hat[j]-f_v_hat[0]*ph_v_hat[j]/ph_v_hat[0];
		  f_hat[j-1]=f_v_hat[j];
		}
		f_v_hat[0]=0;
		normalize(f_hat);
		
		//CONVERT PH_V_HAT AND F_V_HAT FROM FLUID FRAME TO SURF. TETRAD
		//BY JUMPING TO DISK SURFACE AND GOING BACK TO COORD FRAME
		
		if (isbottom == 0) {
		  for (i=0;i<=3;i++) {
		    for (j=0;j<=3;j++) {
		      e_lfs[i][j]=emtop_elf_ik[indexelf(ir_new,jph_new,i,j)];
		    }
		    v_[i]=e_lfs[0][i];
		  }
		  yn[1]=rr[ir_new];
		  yn[2]=emtop_ik[indexr(ir_new,jph_new)];
		  yn[3]=pp[jph_new];
		}
		if (isbottom == 1) {
		  for (i=0;i<=3;i++) {
		    for (j=0;j<=3;j++) {
		      e_lfs[i][j]=embot_elf_ik[indexelf(ir_new,jph_new,i,j)];
		    }
		    v_[i]=e_lfs[0][i];
		  }
		  yn[1]=rr[ir_new];
		  yn[2]=embot_ik[indexr(ir_new,jph_new)];
		  yn[3]=pp[jph_new];
		}
		calc_g(g_dn_ph,g_up_ph,yn);
		po_[0] = g_dn_ph[0][0]*v_[0]+g_dn_ph[0][3]*v_[3];
		po_[1] = g_dn_ph[1][1]*v_[1];
		po_[2] = g_dn_ph[2][2]*v_[2];
		po_[3] = g_dn_ph[0][3]*v_[0]+g_dn_ph[3][3]*v_[3];
		
		calc_e2(e_lf2,w_lf2,g_up_ph,g_dn_ph,po_,v_);
		for (i=0;i<=3;i++) {
		  k_[i]=0;
		  f_[i]=0;
		  for (j=0;j<=3;j++) {
		    k_[i]+=e_lf2[j][i]*ph_v_hat[j];
		    f_[i]+=e_lf2[j][i]*f_v_hat[j];
		  }
		}
		/*
		  printf("\n");
		  printf("%11.4e %11.4e %11.4e %11.4e\n",
		  e_lf2[0][0],e_lf2[0][1],e_lf2[0][2],e_lf2[0][3]);
		  printf("%11.4e %11.4e %11.4e %11.4e\n",
		  e_lf2[1][0],e_lf2[1][1],e_lf2[1][2],e_lf2[1][3]);
		  printf("%11.4e %11.4e %11.4e %11.4e\n",
		  e_lf2[2][0],e_lf2[2][1],e_lf2[2][2],e_lf2[2][3]);
		  printf("%11.4e %11.4e %11.4e %11.4e\n",
		  e_lf2[3][0],e_lf2[3][1],e_lf2[3][2],e_lf2[3][3]);
		*/
		for (i=0;i<=3;i++) {
		  for (j=0;j<=3;j++) {
		    adata1[i+1][j+1]=e_lfs[j][i];
		  }
		}
		ludcmp_js(adata1,4,indx,&dd);
		for (j=0;j<=3;j++) bdata1[j+1]=k_[j];
		lubksb_js(adata1,4,indx,bdata1);
		for (j=0;j<=3;j++) ph_v_hat[j]=bdata1[j+1];
		for (j=0;j<=3;j++) bdata1[j+1]=f_[j];
		lubksb_js(adata1,4,indx,bdata1);
		for (j=0;j<=3;j++) f_v_hat[j]=bdata1[j+1];
		
		E_i = ph_v_hat[0];
		E_i = -yn[4];
		for (j=1;j<=3;j++) {
		  f_v_hat[j]=f_v_hat[j]-f_v_hat[0]*ph_v_hat[j]/ph_v_hat[0];
		  f_hat[j-1]=f_v_hat[j];
		}
		f_v_hat[0]=0;
		normalize(f_hat);	      
		
		//NEW (4/21/09) SCATTERING ALGORITHM ENSURES REFLECTION TO Z>0
		
		//MUST CALCULATE NORMAL TO SURFACE IN LOCAL FRAME
		//THIS SHOULD BE e_z_hat, WHICH CORRESPONDS TO e_lf[3][i]
		
		n_hat[0] = ph_v_hat[1];
		n_hat[1] = ph_v_hat[2];
		n_hat[2] = ph_v_hat[3];
		normalize(n_hat);
		cth0 = n_hat[2];
		//printf("%11.4e %11.4e %11.4e %11.4e\n",
		//     ph_v_hat[0],ph_v_hat[1],ph_v_hat[2],ph_v_hat[3]);
		
		//IF cth0>0 THEN THE PHOTON IS COMING FROM INSIDE THE DISK, AND
		//SHOULD BE LAUNCHED AGAIN WITH RANDOM DIRECTION AND 
		//CHANDRA POLARIZATION
		
		if (cth0 > 0) {
		  //printf("%d %d %d %d %d %d %d %d %10.4e %10.4e %10.4e\n",
		  //     ir,iph,it,ip,ir_old,ir_new,jph_old,jph_new,
		  //     emtop_ik[indexr(ir_new,jph_new)],yn[2],yn[6]);
		  A_fact = 0;
		}
		
		e_x_hat[0]=1;
		e_x_hat[1]=0;
		e_x_hat[2]=0;
		e_y_hat[0]=0;
		e_y_hat[1]=1;
		e_y_hat[2]=0;
		e_z_hat[0]=0;
		e_z_hat[1]=0;
		e_z_hat[2]=1;
		
		//PICK RANDOM SCATTERING ANGLE IN ELF-4 FRAME
		lambda = (double)rand()/(RAND_MAX);
		cth = lambda;
		sth = sqrt(1.-cth*cth);
		lambda = (double)rand()/(RAND_MAX);
		f = lambda*2.*PI;
		
		for (j=0;j<=2;j++) {
		  n_p_hat[j] = cth*e_z_hat[j]+sth*cos(f)*e_x_hat[j]
		    +sth*sin(f)*e_y_hat[j];
		  //NEW PHOTON 4-VELOCITY IN ELF-4 FRAME
		  ph_v_hat[j+1]=n_p_hat[j]*ph_v_hat[0];
		}
		cth = dot(n_p_hat,n_hat);
		mu0 = -dot(n_hat,e_z_hat);
		mu = dot(n_p_hat,e_z_hat);
		ph0 = atan2(n_hat[1],n_hat[0]);
		phi = atan2(n_p_hat[1],n_p_hat[0]);
		
		//NON-ELASTIC SCATTERING- LOSE ENERGY TO ELECTRON
		//(THIS SEEMS NOT QUITE RIGHT FOR DIFFUSE REFLECTION)
		for (iv=0;iv<=Ne;iv++) {
		  rdsh = 1./(1.+nu[iv]*(1./511.)*(1.-cth));
		  //ELASTIC SCATTERING
		  rdsh = 1.;
		  //ESTIMATE FOR ANGLE-AVERAGED 2nd-ORDER SCATTER FOR DIFF. REF.
		  rdsh = 1./(1.+2.*nu[iv]*(1./511.));
		  if (nu[iv]>nu_cut) 
		    Bnu_[iv]=Bnu_[iv]*exp(-pow((nu[iv]-nu_cut)/dnu_cut,2));
		  nu[iv] = nu[iv]*rdsh;
		  dnu[iv] = dnu[iv]*rdsh*rdsh;
		  Bnu_[iv] = Bnu_[iv]/rdsh;
		  //if ((rdsh < 0.2) && (Bnu_[iv] > 1e-10)) 
		  //printf("disk rdsh %g %g %g %g\n",yn[1],nu[iv],Bnu_[iv],rdsh);

		}
		
		//NEW POLARIZATION CALCULATION
		cross(n_hat,e_z_hat,e_parl);
		normalize(e_parl);
		cross(e_parl,n_hat,e_perp);
		
		cross(e_z_hat,n_p_hat,e_parl_f);
		normalize(e_parl_f);
		cross(n_p_hat,e_parl_f,e_perp_f);
		//cross(n_p_hat,e_y_hat,e_parl_f);
		//normalize(e_parl_f);
		//cross(e_parl_f,n_p_hat,e_perp_f);
		cpsi = dot(e_parl,f_hat);
		spsi = dot(e_perp,f_hat);
		I_p = 1;
		Q_p = deg*(2.*cpsi*cpsi-1.);
		I_rr = 0.5*(I_p+Q_p);
		I_ll = 0.5*(I_p-Q_p);
		U_ = deg*2.*spsi*cpsi;
		//printf("%g %g %g %g\n", deg, cpsi, spsi,U_);
		//printf("%g %g %g %g\n", I_ll, I_rr, U_,mu);
		diffuse_reflection(mu0,ph0,mu,phi,&I_ll,&I_rr,&U_);
		I_p = I_ll+I_rr;
		Q_p = I_rr-I_ll;
		U_p = U_;
		//END NEW POL CALC
		
		degp = sqrt(Q_p*Q_p+U_p*U_p)/I_p;
		psip = atan2(U_p,Q_p)/2.;
		for (j=0;j<=2;j++) {
		  //f_v_hat[j+1] = sin(psip)*e_parl_f[j]+cos(psip)*e_perp_f[j];
		  f_v_hat[j+1] = cos(psip)*e_parl_f[j]+sin(psip)*e_perp_f[j];
		}
		E_f = ph_v_hat[0];
		E0 = E_f;

		if (del_eta >= 1) {
		  //PRODUCE IRON FLOURESCENT LINE
		  N_Fe = 0;
		  iv_Fe = -1;
		  E0 = 0.;
		  E_i = 0.;
		  E_f = 0.;
		  E_abs = 9.2;
		  E_em = 7.0;
		  for (iv=0;iv<=Ne;iv++) {
		    if ((nu[iv] > 1.) && (nu[iv] < 7.)) E0+=Bnu_[iv]*dnu[iv]/nu[iv];
		    if ((iv_Fe < 0) && (nu[iv] > E_em)) iv_Fe=iv;
		    if ((nu[iv] > E_abs) && (nu[iv] < 30000.)) {
		      sig_edge = 3.*pow(nu[iv]/E_abs,-3.);
		      sig_edge = sig_edge/(sig_edge+1.);
		      N_Fe += f_Fe*Bnu_[iv]*dnu[iv]/nu[iv]*sig_edge;
		      if ((del_eta == 1)||(del_eta == 2)) 
			Bnu_[iv] = Bnu_[iv]*(1.-sig_edge);
		      //if (Bnu_[iv] > 1e-15) 
		      //printf("%d %g %g %d %g\n",
		      //     iv,nu[iv],Bnu_[iv],iv_Fe,N_Fe*nu[iv_Fe]/dnu[iv_Fe]);
		    }
		    if ((nu[iv] > 1.) && (nu[iv] < 7.)) E_i+=Bnu_[iv]*dnu[iv]/nu[iv];
		    //E_i+=Bnu_[iv]*dnu[iv];
		  }
		  //if (iv_Fe < 0) {
		  //printf("error in iv_Fe %g %g\n",nu[Ne],yn[1]);
		  //iv_Fe = 0;
		  //}
		  if (iv_Fe >= 0) {
		    //printf("made line %d %g %g %g\n",iv_Fe,nu[iv_Fe],yn[1],
		    // (Bnu_[iv_Fe]+N_Fe*nu[iv_Fe]/dnu[iv_Fe])/Bnu_[iv_Fe]);
		    if ((del_eta ==1)||(del_eta ==3))
		      Bnu_[iv_Fe] += N_Fe*nu[iv_Fe]/dnu[iv_Fe];
		    Fe_phot = 1;
		  }
		  for (iv=0;iv<=Ne;iv++) {
		    if ((nu[iv] > 1.) && (nu[iv] < 7.)) E_f+=Bnu_[iv]*dnu[iv]/nu[iv];
		  }
		  //E_f+=Bnu_[iv]*dnu[iv];
		  //if (N_Fe > 1e1) printf("%g %g %g %g %g\n",
		  //		 nu[iv_Fe],E0,E_i,E_f,(E_f-E0)/E0);
		  //nu[iv_Fe],Bnu_[iv_Fe-1],Bnu_[iv_Fe],Bnu_[iv_Fe+1],
		  //   N_Fe*nu[iv_Fe]/dnu[iv_Fe]);
		  //END IRON LINE 
		}
		//CONVERT BACK FROM DISK FRAME TO COORDINATE BASIS
		for (i=0;i<=3;i++) {
		  ph_v_p[i]=0;
		  f_[i]=0;
		  for (j=0;j<=3;j++) {
		    ph_v_p[i]+=e_lfs[j][i]*ph_v_hat[j];
		    f_[i]+=e_lfs[j][i]*f_v_hat[j];
		  }
		}
		deg = degp;
		A_fact = A_fact*I_p;
		kap1i = aa*cos(t)*((ph_v_p[0]*f_[1]-ph_v_p[1]*f_[0])
				   +aa*sin(t)*sin(t)*(ph_v_p[1]*f_[3]-ph_v_p[3]*f_[1]))
		  +r*sin(t)*((r*r+aa*aa)*(ph_v_p[3]*f_[2]-ph_v_p[2]*f_[3])
			     -aa*(ph_v_p[0]*f_[2]-ph_v_p[2]*f_[0]));
		kap2i = r*(ph_v_p[0]*f_[1]-ph_v_p[1]*f_[0]
			   +aa*sin(t)*sin(t)*(ph_v_p[1]*f_[3]-ph_v_p[3]*f_[1]))-
		  aa*sin(t)*cos(t)*((r*r+aa*aa)*(ph_v_p[3]*f_[2]-ph_v_p[2]*f_[3])
				    -aa*(ph_v_p[0]*f_[2]-ph_v_p[2]*f_[0]));
		
		yn[4] = g_dn_ph[0][0]*ph_v_p[0]+g_dn_ph[0][3]*ph_v_p[3];
		yn[5] = g_dn_ph[1][1]*ph_v_p[1];
		yn[6] = g_dn_ph[2][2]*ph_v_p[2];
		yn[7] = g_dn_ph[0][3]*ph_v_p[0]+g_dn_ph[3][3]*ph_v_p[3];
		E_f = -yn[4];
		pdv_em = (ph_v_p[0]*po_[0]+ph_v_p[1]*po_[1]+
			  ph_v_p[2]*po_[2]+ph_v_p[3]*po_[3]);
		
		if ((E_f/E_i > 5e5)&&(E_f/E_i < 5e5)) {
		  printf("Ef/Ei = %12.5e %12.5e %12.5e %12.5e %12.5e\n",
			 E_f/E_i, E_i,yn[3],r,cth);
		}
		
		calc_e4(e_lfs,w_lfs,g_up_ph,g_dn_ph,
			po_,v_,dx_r,dx_p,isbottom,&G_fact);
		if (A_fact > 0) {
		  Ndisk++;
		  dscat++;
		}
		just_scattered=1; 
		//TO COMPARE WITH SPECTRUM.C, MAKE MIDPLANE OPAQUE
		//A_fact = 0;

		//re-normalize to ensure p_.p_=0
		calc_g(g_dn_ph,g_up_ph,yn);
		za = g_up_ph[0][0];
		zb = 2.*yn[7]*g_up_ph[0][3];
		zc = yn[5]*yn[5]*g_up_ph[1][1]+
		  yn[6]*yn[6]*g_up_ph[2][2]+yn[7]*yn[7]*g_up_ph[3][3];
		zq = zb*zb-4*za*zc;
		if (zq < 0) printf("normalization error\n");
		if (zq > 0) {
		  s1 = (-zb+sqrt(zq))/(2.*za);
		  s2 = (-zb-sqrt(zq))/(2.*za);
		  if (fabs(s1-yn[4]) < fabs(s2-yn[4])) yn[4]= s1;
		  if (fabs(s1-yn[4]) > fabs(s2-yn[4])) yn[4]= s2;
		}
		for (j=0;j<=3;j++) p_[j]=yn[j+4];
		E_f = dot_g4(g_up_ph,p_,p_);
		//if (fabs(E_f) > 1e-20) printf("a %g\n",E_f);
		for (j=0;j<=7;j++) {
		  y2[j]=yn[j];
		}
	      }	  /*********INTERSECTION WITH DISK**********/
	    }

	    //update the photon position+momentum 
	    for (j=0;j<=7;j++) {
	      yn[j]=y2[j];
	    }	

	    /****************PASSING THROUGH Z-AXIS********/
	    if (yn[2] < 0) {
	      yn[2] = -yn[2];
	      yn[3] = yn[3]+PI;
	      yn[6] = -yn[6];
	    }
	    if (yn[2] > PI) {
	      yn[2] = 2.*PI-yn[2];
	      yn[3] = yn[3]+PI;
	      yn[6] = -yn[6];
	    }
	    	    
	    //ALWAYS UPDATE THE COORDINATE TIME
	    yn[0]=y[0]+dt;
	    
	    /*****************CORONAL SCATTERING*****************/
	    if ((just_scattered == 0)&&(yn[1] > 1.02*Rhor)&&(tau_es > 10)) {
	      //calc_g(g_dn_ph,g_up_ph,yn);
	      r = yn[1];
	      r2 = r*r;
	      cth2 = cos(yn[2])*cos(yn[2]);
	      sth2 = 1.-cth2;
	      calc_g(g_dn_ph,g_up_ph,yn);
	      //g_dn_ph[2][2]=r2+a2*cth2;
	      //g_dn_ph[1][1]=g_dn_ph[2][2]/(r2-2*M*r+a2);
	      //g_dn_ph[3][3]=(r2+a2+2*M*a2*r*sth2/g_dn_ph[2][2])*sth2;
	      dtau = sqrt(g_dn_ph[1][1]*(yn[1]-y[1])*(yn[1]-y[1])
			  +g_dn_ph[2][2]*(yn[2]-y[2])*(yn[2]-y[2])
			  +g_dn_ph[3][3]*(yn[3]-y[3])*(yn[3]-y[3]));
	      ph_v[0] = g_up_ph[0][0]*yn[4]+g_up_ph[0][3]*yn[7];
	      ph_v[1] = g_up_ph[1][1]*yn[5];
	      ph_v[2] = g_up_ph[2][2]*yn[6];
	      ph_v[3] = g_up_ph[0][3]*yn[4]+g_up_ph[3][3]*yn[7];
	      for (j=0;j<=3;j++) {
		part_x[j]=yn[j];
		k_[j]=ph_v[j];
		p_[j]=yn[j+4];
	      }

	      //check to see if re-normalization necessary
	      E_f = dot_g4(g_up_ph,p_,p_);
	      if (fabs(E_f) > 1e-10) {
		za = g_up_ph[0][0];
		zb = 2.*yn[7]*g_up_ph[0][3];
		zc = yn[5]*yn[5]*g_up_ph[1][1]+
		  yn[6]*yn[6]*g_up_ph[2][2]+yn[7]*yn[7]*g_up_ph[3][3];
		zq = zb*zb-4*za*zc;
		if (zq < 0) printf("normalization error\n");
		if (zq > 0) {
		  s1 = (-zb+sqrt(zq))/(2.*za);
		  s2 = (-zb-sqrt(zq))/(2.*za);
		  if (fabs(s1-yn[4]) < fabs(s2-yn[4])) yn[4]= s1;
		  if (fabs(s1-yn[4]) > fabs(s2-yn[4])) yn[4]= s2;
		}
	      }
	      
	      lookup_data(part_x,rr,tt,pp,g_dn_ph,
			  rho_ijk,T_ijk,bb_ijk,ut_ijk,ur_ijk,uz_ijk,up_ijk,
			  weights,&rho,&T_e,&Bmag,v_);
	      E0 = g_dn_ph[0][0]*v_[0]*v_[0]+2.*g_dn_ph[0][3]*v_[0]*v_[3]
		+g_dn_ph[1][1]*v_[1]*v_[1]+g_dn_ph[2][2]*v_[2]*v_[2]
		+g_dn_ph[3][3]*v_[3]*v_[3];
	      
	      rho_bar = rho;
	      dtau_es = kapp_es*rho*dtau*R_g;
	      tau_tot=tau_tot+dtau_es;
	      l_tot = l_tot+dtau*R_g;
	      lambda = (double)rand()/(RAND_MAX);
	      if (E0 > -0.99) printf("vcor %ld %ld %g %g %g %g %g\n",
				     steps,iscat,dtau,rho,r,yn[2],dtau_es);

	      //cyclotron self-absorption
	      if (ibottom >= 20) {
		//rdsh is redshift of photon packet in corona rest frame,
		//relative to the point of last measurement
		for (j=0;j<=3;j++) p_[j]=yn[j+4];
		rdsh = (v_[0]*p_[0]+v_[1]*p_[1]+v_[2]*p_[2]+v_[3]*p_[3])/
		  pdv_em;
		for (iv=0;iv<=Ne;iv++) {
		  nu[iv] = nu[iv]*rdsh;
		  dnu[iv] = dnu[iv]*rdsh;
		}
		pdv_em = (v_[0]*p_[0]+v_[1]*p_[1]+v_[2]*p_[2]+v_[3]*p_[3]);
		
		gam2 = pow((kB_ev*T_e/5.11e5+1.),2.);
		nuc = qe*sqrt(Bmag)/(2*PI*me*cc); //nu_c in Hz
		nuc = nuc/2.42e17; //nu_c in keV
		if (nuc == 0) nuc = 1.e-5;
		n_e = rho/mp;
		P0 = lookup_Pnorm(Pnorm,nu0,T_e*kB_ev*1e-3,Ne);
		//printf("%g %g\n",T_e,P0);
		if (P0 > 0) P0 = 4./3.*sigma_T*cc*(gam2-1)*(Bmag/(8.*PI))*n_e/
		  (nuc*2.42e17*P0*4.*PI);
		for (iv=0;iv<=Ne;iv++) {
		  //blackbody brightness in units of [ergs/cm^2/s/Hz/Sr]	      
		  B_nu[iv] = 0.;
		  if (T_e > 0) {
		    x = nu[iv]/(kB_ev*T_e*1e-3);
		    B_nu[iv]=2.08e5*(nu[iv]*nu[iv]*nu[iv])/(exp(x)-1.);
		  }
		  //synchrotron emissivity in units of [ergs/cm^3/s/Hz/Sr]
		  j_nu[iv]= 0.;
		  x = nu[iv]/nuc;
		  if ((T_e > 1e7)&&(Bmag > 0)&&(x <= 1)) 
		    j_nu[iv]=P0*F_cycl(x,kB_ev*T_e*1e-3);
		  //synchrotron self-absorption in units of [1/cm]
		  alpha_nu[iv] = 0;
		  if (B_nu[iv] > 0) alpha_nu[iv]=j_nu[iv]/B_nu[iv];
		  //Brightness absorbed in the corona cell
		  B_nu[iv]=Bnu_[iv]*(1.-exp(-alpha_nu[iv]*dtau*R_g)); 
		  //Brightness exiting the corona cell
		  //Bnu_[iv]=Bnu_[iv]*exp(-alpha_nu[iv]*dtau*R_g);
		  Bnu_[iv]=Bnu_[iv]-B_nu[iv];
		  //printf("%g %g %g %g %g\n",
		  //       nu[iv],dtau,alpha_nu[iv],j_nu[iv],B_nu[iv]);
		  //power absorbed by the cell is negative, since it looks like 
		  //negative cooling, i.e., heating of the corona
		  /*
		    for (i=0;i<=1;i++) {
		    for (j=0;j<=1;j++) {
		    for (k=0;k<=1;k++) {
		    wght = weights[i+6]*weights[j+8]*weights[k+10];
		    corpow_ijk[indexijk((int)weights[i],(int)weights[j+2],
		    (int)weights[k+4])] 
		    -= wght*B_nu[iv]*dnu[iv]*A_fact;
		    }
		    }
		    }
		  */
		  corpow_ijk[indexijk((int)weights[0],(int)weights[2],
				      (int)weights[4])] 
		    -= wght*B_nu[iv]*dnu[iv]*A_fact;
		  if (nu[iv] > nuc) {
		    iv = Ne+1;
		  }
		}
	      }
	      //synchrotron self-absorption
	      if (ibottom >= 20) {
		//rdsh is redshift of photon packet in corona rest frame,
		//relative to the point of last measurement
		for (j=0;j<=3;j++) p_[j]=yn[j+4];
		rdsh = (v_[0]*p_[0]+v_[1]*p_[1]+v_[2]*p_[2]+v_[3]*p_[3])/
		  pdv_em;
		for (iv=0;iv<=Ne;iv++) {
		  nu[iv] = nu[iv]*rdsh;
		  dnu[iv] = dnu[iv]*rdsh;
		}
		pdv_em = (v_[0]*p_[0]+v_[1]*p_[1]+v_[2]*p_[2]+v_[3]*p_[3]);
		
		gam2 = pow((kB_ev*T_e/5.11e5+1.),2.);
		wc = 3.*PI/8.*gam2*qe*sqrt(Bmag)/(me*cc); //omega_c in Hz
		wc = wc/2.42e17; //omega_c in keV
		n_e = rho/mp;
		P0 = 4./sqrt(3.)/(PI*PI)*(qe*qe*qe/(me*cc*cc))* 
		  sqrt(Bmag)*(1.-1./gam2)*n_e/(4.*PI);
		if (wc == 0) wc = 1.e-5;
		for (iv=0;iv<=Ne;iv++) {
		  //blackbody brightness in units of [ergs/cm^2/s/Hz/Sr]	      
		  B_nu[iv] = 0.;
		  if (T_e > 0) {
		    x = nu[iv]/(kB_ev*T_e*1e-3);
		    B_nu[iv]=2.08e5*(nu[iv]*nu[iv]*nu[iv])/(exp(x)-1.);
		  }
		  //synchrotron emissivity in units of [ergs/cm^3/s/Hz/Sr]
		  j_nu[iv]= 0.;
		  x = nu[iv]/wc;
		  if ((T_e > 0)&&(Bmag > 0)) 
		    j_nu[iv]=P0*F_sync(x);
		  //synchrotron self-absorption in units of [1/cm]
		  alpha_nu[iv] = 0;
		  if (B_nu[iv] > 0) alpha_nu[iv]=j_nu[iv]/B_nu[iv];
		  //Brightness absorbed in the corona cell
		  B_nu[iv]=Bnu_[iv]*(1.-exp(-alpha_nu[iv]*dtau*R_g)); 
		  //Brightness exiting the corona cell
		  //Bnu_[iv]=Bnu_[iv]*exp(-alpha_nu[iv]*dtau*R_g);
		  Bnu_[iv]=Bnu_[iv]-B_nu[iv];
		  //printf("%g %g %g %g %g\n",
		  //       nu[iv],dtau,alpha_nu[iv],j_nu[iv],B_nu[iv]);
		  //power absorbed by the cell is negative, since it looks like 
		  //negative cooling, i.e., heating of the corona
		  /*
		    for (i=0;i<=1;i++) {
		    for (j=0;j<=1;j++) {
		    for (k=0;k<=1;k++) {
		    wght = weights[i+6]*weights[j+8]*weights[k+10];
		    corpow_ijk[indexijk((int)weights[i],(int)weights[j+2],
		    (int)weights[k+4])] 
		    -= wght*B_nu[iv]*dnu[iv]*A_fact;
		    }
		    }
		    }
		  */
		  corpow_ijk[indexijk((int)weights[0],(int)weights[2],
				      (int)weights[4])] 
		    -= wght*B_nu[iv]*dnu[iv]*A_fact;
		  if (alpha_nu[iv]*dtau*R_g < 1e-4) {
		    iv = Ne+1;
		  }
		}
	      }
	      //if ((it == 14)&&(ip == 0)) {
	      //calc_g(g_dn_ph,g_up_ph,yn);
	      //printf("%g %g %g\n",rdsh,pdv_em,dot_g4(g_up_ph,p_,p_));
	      //}
	      //printf("%d %d %g\n",it, steps, lambda);
	      if (lambda < (1.-exp(-dtau_es))) {
		//printf("%10.4e %10.4e %10.4e %10.4e\n",r,yn[2],yn[3],dtau);
		//printf("%d %g %g %g %g\n",steps,rho,dtau_es,T_e,R_g);
		pow_i = 0;
		rdsh0 = p_[0]/pdv_em;
		for (iv=0;iv<=Ne;iv++) pow_i+=Bnu_[iv]*dnu[iv];
		//pow_i = pow_i*rdsh0;
		calc_g(g_dn_ph,g_up_ph,yn);
		iscat++;
		cth = cos(yn[2]);
		sth = sin(yn[2]);
		r2=r*r;
		ph_v[0] = g_up_ph[0][0]*yn[4]+g_up_ph[0][3]*yn[7];
		ph_v[1] = g_up_ph[1][1]*yn[5];
		ph_v[2] = g_up_ph[2][2]*yn[6];
		ph_v[3] = g_up_ph[0][3]*yn[4]+g_up_ph[3][3]*yn[7];
		
		for (j=0;j<=3;j++) {
		  p_[j]=yn[j+4];
		  po_[j]=0;
		  k_[j]=ph_v[j];
		}
		po_[0] = g_dn_ph[0][0]*v_[0]+g_dn_ph[0][3]*v_[3];
		po_[1] = g_dn_ph[1][1]*v_[1];
		po_[2] = g_dn_ph[2][2]*v_[2];
		po_[3] = g_dn_ph[0][3]*v_[0]+g_dn_ph[3][3]*v_[3];

		//rdsh is redshift of scattered photons in corona rest frame
		rdsh = (v_[0]*p_[0]+v_[1]*p_[1]+v_[2]*p_[2]+v_[3]*p_[3])/
		  pdv_em;
		for (iv=0;iv<=Ne;iv++) {
		  nu[iv] = nu[iv]*rdsh;
		  dnu[iv] = dnu[iv]*rdsh;
		}
		
		bdata1[1]=kap1i;
		bdata1[2]=kap2i;
		bdata1[3]=0;
		bdata1[4]=0;
		adata1[1][1]=-aa*cth*k_[1]+r*sth*aa*k_[2];
		adata1[1][2]=aa*cth*k_[0]-aa*sth*sth*k_[3];
		adata1[1][3]=r*sth*((r2+a2)*k_[3]-aa*k_[0]);
		adata1[1][4]=a2*cth*sth*sth*k_[1]+r*sth*(-(r2+a2)*k_[2]);
		adata1[2][1]=r*(-k_[1])+a2*sth*cth*k_[2];
		adata1[2][2]=r*(k_[0]-aa*sth*sth*k_[3]);
		adata1[2][3]=aa*sth*cth*((r2+a2)*k_[3]-aa*k_[0]);
		adata1[2][4]=r*aa*sth*sth*k_[1]-aa*sth*cth*(-(r2+a2)*k_[2]);
		for (j=0;j<=3;j++) {
		  adata1[3][j+1]=p_[j];
		  adata1[4][j+1]=po_[j];
		}
		
		ludcmp_js(adata1,4,indx,&dd);
		lubksb_js(adata1,4,indx,bdata1);
		for (j=0;j<=3;j++) f_[j]=bdata1[j+1];
		
		//transform into local gas frame
		//calc_e(e_lfs,w_lfs,g_up_ph,g_dn_ph,po_,v_);
		calc_e2(e_lfs,w_lfs,g_up_ph,g_dn_ph,po_,v_);
		for (i=0;i<=3;i++) {
		  for (j=0;j<=3;j++) {
		    adata1[i+1][j+1]=e_lfs[j][i];
		  }
		}
		E0 = g_dn_ph[0][0]*k_[0]*k_[0]+2.*g_dn_ph[0][3]*k_[0]*k_[3]
		  +g_dn_ph[1][1]*k_[1]*k_[1]+g_dn_ph[2][2]*k_[2]*k_[2]
		  +g_dn_ph[3][3]*k_[3]*k_[3];
		//printf("a %d %g\n",steps,E0);
		ludcmp_js(adata1,4,indx,&dd);
		for (j=0;j<=3;j++) bdata1[j+1]=k_[j];
		lubksb_js(adata1,4,indx,bdata1);
		for (j=0;j<=3;j++) ph_v_hat[j]=bdata1[j+1];
		for (j=0;j<=3;j++) bdata1[j+1]=f_[j];
		lubksb_js(adata1,4,indx,bdata1);
		for (j=0;j<=3;j++) f_v_hat[j]=bdata1[j+1];
		E0 = -ph_v_hat[0]*ph_v_hat[0]+ph_v_hat[1]*ph_v_hat[1]
		  +ph_v_hat[2]*ph_v_hat[2]+ph_v_hat[3]*ph_v_hat[3];
		//printf("ap %g\n",E0);

		//printf("pre-scatter p0  %12.5g\n",ph_v_hat[0]);
		E_i = ph_v_hat[0];
		//PICK RANDOM 3-VELOCITY FOR ELECTRON IN CORONA FRAME
		lambda = (double)rand()/(RAND_MAX);
		//ELECTRON TEMP IN MeV
		T_e = T_e/(11000.)/1.0e6;
		//BULK COMPTONIZATION
		if (T_cor == 0) T_e = 1e-6;
		Z1 = -2.*T_e/0.511*log(1.-lambda);
		Z2 = Z1*Z1;
		beta = sqrt(0.5*(-Z2+sqrt(Z2*Z2+4.*Z2)));

		cth = 2.0*(double)rand()/(RAND_MAX)-1.0;
		sth = sqrt(1.-cth*cth);
		lambda = (double)rand()/(RAND_MAX);
		f = lambda*2.*PI;
		n_hat[0] = sth*cos(f);
		n_hat[1] = sth*sin(f);
		n_hat[2] = cth;
		boost(beta,n_hat,ph_v_hat);
		boost(beta,n_hat,f_v_hat);
		
		rdsh = ph_v_hat[0]/E_i;
		for (j=0;j<=3;j++) ph_v_hat[j]=ph_v_hat[j]/rdsh;
		for (iv=0;iv<=Ne;iv++) {
		  nu[iv] = nu[iv]*rdsh;
		  dnu[iv] = dnu[iv]*rdsh;
		}

		//RE-WRITE f_ SO THAT f_[0]=0, MAINTAINING f_.k_=0 
		//(CONNERS, PIRAN, AND STARK 1980) 
		
		for (j=1;j<=3;j++) {
		  f_v_hat[j]=f_v_hat[j]-f_v_hat[0]*ph_v_hat[j]/ph_v_hat[0];
		  f_hat[j-1]=f_v_hat[j];
		}
		f_v_hat[0]=0;
		normalize(f_hat);
		//printf("post %g %g %g %g\n", 
		//   ph_v_hat[0],ph_v_hat[1],ph_v_hat[2],ph_v_hat[3]);
		
		//CALCULATE THOMSON SCATTERING ANGLE IN ELECTRON FRAME
		e_z_hat[0] = ph_v_hat[1];
		e_z_hat[1] = ph_v_hat[2];
		e_z_hat[2] = ph_v_hat[3];
		normalize(e_z_hat);
		r_hat[0] = (double)rand()/(RAND_MAX)-0.5;
		r_hat[1] = (double)rand()/(RAND_MAX)-0.5;
		r_hat[2] = (double)rand()/(RAND_MAX)-0.5;
		cross(e_z_hat,r_hat,e_x_hat);
		normalize(e_x_hat);
		cross(e_z_hat,e_x_hat,e_y_hat);
		
		/*OLD scattering algorithm selects theta and psi from their
		  actual probability distribution */
	      
		calc_scat_angles(deg,&cth,&f);
		sth = sqrt(1.-cth*cth);
		cpsi = cos(f);
		spsi = sin(f);
		//NON-ELASTIC SCATTERING- LOSE ENERGY TO ELECTRON
		for (iv=0;iv<=Ne;iv++) {
		  rdsh = 1./(1.+nu[iv]*(1./511.)*(1.-cth));
		  if (nu[iv]>nu_cut) 
		    Bnu_[iv]=Bnu_[iv]*exp(-pow((nu[iv]-nu_cut)/dnu_cut,2));
		  nu[iv] = nu[iv]*rdsh;
		  dnu[iv] = dnu[iv]*rdsh*rdsh;
		  Bnu_[iv] = Bnu_[iv]/rdsh;
		  //if ((rdsh < 0.2) && (Bnu_[iv] > 1e-10)) 
		  //printf("cor rdsh %g %g %g %g\n",yn[1],nu[iv],Bnu_[iv],rdsh);
		}
		
		cross(e_z_hat,f_hat,fp_hat);
		normalize(fp_hat);
		for (j=0;j<=2;j++) {
		  e_x_hat[j]=cpsi*f_hat[j]+spsi*fp_hat[j];
		}
		normalize(e_x_hat);
		cross(e_z_hat,e_x_hat,e_y_hat);
		for (j=0;j<=2;j++) {
		  n_p_hat[j] = cth*e_z_hat[j]+sth*e_x_hat[j];
		  ph_v_hat[j+1]=n_p_hat[j]*ph_v_hat[0];
		}
		E0 = -ph_v_hat[0]*ph_v_hat[0]+ph_v_hat[1]*ph_v_hat[1]
		  +ph_v_hat[2]*ph_v_hat[2]+ph_v_hat[3]*ph_v_hat[3];
		//E_f = n_p_hat[0]*n_p_hat[0]+n_p_hat[1]*n_p_hat[1]+n_p_hat[2]*n_p_hat[2];
		//E_f = f_hat[0]*e_z_hat[0]+f_hat[1]*e_z_hat[1]+f_hat[2]*e_z_hat[2];
		if (fabs(E0) > 1e-2) printf("b %g %g %g\n",E0, E_f, E_i);
		cross(e_z_hat,n_p_hat,e_perp);
		normalize(e_perp);
		cross(e_perp,e_z_hat,e_parl);
		for (j=0;j<=2;j++) e_perp_f[j]=e_perp[j];
		cross(e_perp_f,n_p_hat,e_parl_f);
		cpsi = dot(e_parl,f_hat);
		spsi = dot(e_perp,f_hat);
		
		I_p = 1;
		Q_p = deg*(2.*cpsi*cpsi-1.);
		I_rr = 0.5*(I_p+Q_p);
		I_ll = 0.5*(I_p-Q_p);
		U_ = deg*2.*spsi*cpsi;
		
		I_p = 3./2.*(cth*cth*I_rr+I_ll);
		Q_p = 3./2.*(cth*cth*I_rr-I_ll);
		U_p = 3./2.*cth*U_;
		degp = sqrt(Q_p*Q_p+U_p*U_p)/I_p;
		psip = atan2(U_p,Q_p)/2.;
		for (j=0;j<=2;j++) {
		  f_v_hat[j+1] = cos(psip)*e_parl_f[j]+sin(psip)*e_perp_f[j];
		}
		
		//BOOST BACK INTO CORONA FRAME, THEN CONVERT TO COORDINATE BASIS
		beta = -beta;
		//printf("post-scatter p0  %12.5g\n",ph_v_hat[0]);
		E_f = ph_v_hat[0];
		E0 = -ph_v_hat[0]*ph_v_hat[0]+ph_v_hat[1]*ph_v_hat[1]
		  +ph_v_hat[2]*ph_v_hat[2]+ph_v_hat[3]*ph_v_hat[3];
		if (fabs(E0) > 1e-4) {
		  printf("neg energy %ld %ld %ld %ld %g %g %g %g %g %g\n",it,ip,ir,iph,
			 E0,dot_g4(g_up_ph,p_,p_),dot_g4(g_dn_ph,v_,v_),yn[1],yn[2],rdsh);
		  sleep(1);
		}
		boost(beta,n_hat,ph_v_hat);
		boost(beta,n_hat,f_v_hat);
		rdsh = ph_v_hat[0]/E_f;

		if (rdsh < 0) {
		  printf("%g %g %g %g %g %g %g %g\n",rdsh,ph_v_hat[0],E_f,
			 //beta,n_hat[0],n_hat[1],n_hat[2]);
			 v_[0],v_[1],v_[2],v_[3],dot_g4(g_dn_ph,v_,v_));
		  E_f = -ph_v_hat[0]*ph_v_hat[0]+ph_v_hat[1]*ph_v_hat[1]
		    +ph_v_hat[2]*ph_v_hat[2]+ph_v_hat[3]*ph_v_hat[3];
		  //printf("%g %g\n",E0,E_f);
		}
		for (j=0;j<=3;j++) ph_v_hat[j]=ph_v_hat[j]/rdsh;
		for (iv=0;iv<=Ne;iv++) {
		  nu[iv] = nu[iv]*rdsh;
		  dnu[iv] = dnu[iv]*rdsh;
		}

		for (j=1;j<=3;j++) {
		  f_v_hat[j]=f_v_hat[j]-f_v_hat[0]*ph_v_hat[j]/ph_v_hat[0];
		}
		f_v_hat[0]=0;
		
		E_f = ph_v_hat[0];
		E0 = E_f;
		//printf("Ef/Ei = %12.5g %12.5g %12.5g\n",E_f/E_i, E_f,A_fact);
		//fprintf(outfile,"%12.5g\n",Z1);
		
		for (i=0;i<=3;i++) {
		  ph_v_p[i]=0;
		  f_[i]=0;
		  for (j=0;j<=3;j++) {
		    ph_v_p[i]+=e_lfs[j][i]*ph_v_hat[j];
		    f_[i]+=e_lfs[j][i]*f_v_hat[j];
		  }
		}
		deg = degp;
		//weight photon by angle-dependent differential cross-section
		//A_fact = A_fact*I_p;
		kap1i = aa*cos(t)*((ph_v_p[0]*f_[1]-ph_v_p[1]*f_[0])
				   +aa*sin(t)*sin(t)*(ph_v_p[1]*f_[3]-ph_v_p[3]*f_[1]))
		  +r*sin(t)*((r*r+aa*aa)*(ph_v_p[3]*f_[2]-ph_v_p[2]*f_[3])
			     -aa*(ph_v_p[0]*f_[2]-ph_v_p[2]*f_[0]));
		kap2i = r*(ph_v_p[0]*f_[1]-ph_v_p[1]*f_[0]
			   +aa*sin(t)*sin(t)*(ph_v_p[1]*f_[3]-ph_v_p[3]*f_[1]))-
		  aa*sin(t)*cos(t)*((r*r+aa*aa)*(ph_v_p[3]*f_[2]-ph_v_p[2]*f_[3])
				    -aa*(ph_v_p[0]*f_[2]-ph_v_p[2]*f_[0]));
		
		yn[4] = g_dn_ph[0][0]*ph_v_p[0]+g_dn_ph[0][3]*ph_v_p[3];
		yn[5] = g_dn_ph[1][1]*ph_v_p[1];
		yn[6] = g_dn_ph[2][2]*ph_v_p[2];
		yn[7] = g_dn_ph[0][3]*ph_v_p[0]+g_dn_ph[3][3]*ph_v_p[3];
		pdv_em = (ph_v_p[0]*po_[0]+ph_v_p[1]*po_[1]+
			  ph_v_p[2]*po_[2]+ph_v_p[3]*po_[3]);
		pow_f = 0;
		for (iv=0;iv<=Ne;iv++) pow_f+=Bnu_[iv]*dnu[iv];
		for (i=0;i<=1;i++) {
		  for (j=0;j<=1;j++) {
		    for (k=0;k<=1;k++) {
		      wght = weights[i+6]*weights[j+8]*weights[k+10];
		      corpow_ijk[indexijk((int)weights[i],(int)weights[j+2],
					  (int)weights[k+4])] 
			+= wght*(pow_f-pow_i)*A_fact;
		    }
		  }
		}
		//re-normalize to ensure p_.p_=0
		//calc_g(g_dn_ph,g_up_ph,yn);
		for (j=0;j<=3;j++) p_[j]=yn[j+4];
		E_f = dot_g4(g_up_ph,p_,p_);

		if (fabs(E_f) > 1e-10) {
		  za = g_up_ph[0][0];
		  zb = 2.*yn[7]*g_up_ph[0][3];
		  zc = yn[5]*yn[5]*g_up_ph[1][1]+
		    yn[6]*yn[6]*g_up_ph[2][2]+yn[7]*yn[7]*g_up_ph[3][3];
		  zq = zb*zb-4*za*zc;
		  if (zq < 0) printf("normalization error\n");
		  if (zq > 0) {
		    s1 = (-zb+sqrt(zq))/(2.*za);
		    s2 = (-zb-sqrt(zq))/(2.*za);
		    if (fabs(s1-yn[4]) < fabs(s2-yn[4])) yn[4]= s1;
		    if (fabs(s1-yn[4]) > fabs(s2-yn[4])) yn[4]= s2;
		  }
		}
		Fe_phot = 0;
	      }
	    }
	    /*************END CORONAL SCATTERING*****************/
	    
	    if (yn[1] < 1.02*Rhor) {
	      wlo = (yn[1]-1.02*Rhor)/(yn[1]-y[1]);
	      whi = (1.02*Rhor-y[1])/(yn[1]-y[1]);
	      for (j=0;j<=7;j++) {
		yn[j] = wlo*y[j]+whi*yn[j];
	      }
	      rdsh = yn[4]/pdv_em;
	      //printf("%12.5e %12.5e %12.5e\n",rr[ir],yn[4],pdv_em);
	      for (iv=0;iv<=Ne;iv++) {
		nu[iv] = nu0[iv]*rdsh;
		dnu[iv] = dnu0[iv]*rdsh;
		//linear scale
		if (spec_model == 1) 
		  je=(int)floor((nu[iv]-e_min)/(e_max-e_min)*Ne);
		//log scale
		if (spec_model == 2) 
		  je=(int)floor((Ne)*log(nu[iv]/e_min)/log(e_max/e_min));
		if ((je >= 0)&(je < Ne)) {
		  elo = je;
		  ehi = je+1;
		  Eph = nu[iv];
		  wlo = (nu0[ehi]-Eph)/(nu0[ehi]-nu0[elo]);
		  whi = (Eph-nu0[elo])/(nu0[ehi]-nu0[elo]);
		  if (wlo < 0) {
		    printf("%ld %ld %g %g %g\n",iv,je,nu[iv],dnu[iv],Bnu_[iv]);
		  }
		  Cspecr[indexre(ir,elo)] += wlo*A_fact*Bnu_[iv];
		  Cspecr[indexre(ir,ehi)] += whi*A_fact*Bnu_[iv];
		  if (iv == -50) {
		    printf("%ld %ld %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
			   iv,je,nu0[elo],Eph,nu0[ehi],p_[0],rdsh,wlo+whi);
		  }
		}
	      }	    
	      yn[1]=1.01*Rhor;
	      Ncapt += 1.;
	      //printf("%g %g %d %d %12.5e\n", Nescape, Ncapt,it, ip, yn[1]);	  
	    }
	    if (yn[1] > Rshell) {
	      //printf("%d %d %d %d\n", it,ip,iscat,dscat);
	      //if ((ir == 10)&&(it == 126)&&(ip == 124)) printf("OK");
	      wlo = (yn[1]-Rshell)/(yn[1]-y[1]);
	      whi = (Rshell-y[1])/(yn[1]-y[1]);
	      for (j=0;j<=7;j++) {
		yn[j] = wlo*y[j]+whi*yn[j];
	      }
	      jth = (int)(fabs(cos(yn[2]))*(Nth_obs+1.));
	      if (TWO_SIDED >= 1) {
		jth = (int)((cos(yn[2])+1.)/2.*(Nth_obs+1.));
	      }
	      if (jth > Nth_obs) jth = Nth_obs;
		
	      while (yn[3] < 0) yn[3]+=2.*PI;
	      jph = (int)(fmod(yn[3],2.*PI)/(2.*PI)*(Nph_obs+1));
	      if (jph > Nph_obs) jph = Nph_obs;
	      jt = (int)(fabs(yn[0]-Tmin)/(Tmax-Tmin)*Nt);
	      if (jt > Nt) jt=Nt;
	      Iobs[indexph(jth,jph,jt)]=0;
	      if (imageprint == 1) {
		e_x[0]=-sin(yn[3]);
		e_x[1]=cos(yn[3]);
		e_x[2]=0;
		e_y[0]=-cos(yn[2])*cos(yn[3]);
		e_y[1]=-cos(yn[2])*sin(yn[3]);
		e_y[2]=sin(yn[2]);
		x_[0]=yn[1];
		x_[1]=yn[2];
		x_[2]=yn[3];
		p_hat[0]=yn[5];
		p_hat[1]=yn[6]/(yn[1]*yn[1]);
		p_hat[2]=yn[7]/(yn[1]*yn[1]*(sin(yn[2])*sin(yn[2])+1e-5));
		spherv_to_cartv(x_,p_hat,r_hat);
		normalize(r_hat);
		ix = (int)((dot(r_hat,e_x)/(fov/Rshell)+1)*Ni/2.);
		iy = (int)((dot(r_hat,e_y)/(fov/Rshell)+1)*Ni/2.);
	      }

	      //Polarization
	      calc_g(g_dn_ph,g_up_ph,yn);
	      r = yn[1];
	      t = yn[2];
	      r2 = r*r;
	      k_[0]=g_up_ph[0][0]*yn[4]+g_up_ph[0][3]*yn[7];
	      k_[1]=g_up_ph[1][1]*yn[5];
	      k_[2]=g_up_ph[2][2]*yn[6];
	      k_[3]=g_up_ph[3][0]*yn[4]+g_up_ph[3][3]*yn[7];
	      for (j=0;j<=3;j++) {
		p_[j]=yn[j+4];
		po_[j]=0;
		v_[j]=0;
	      }
	      //momentum, velocity of ZAMO at Rshell
	      po_[0]=-1/sqrt(-g_up_ph[0][0]);
	      v_[0]=g_up_ph[0][0]*po_[0];
	      v_[3]=g_up_ph[0][3]*po_[0];
	      bdata1[1]=kap1i;
	      bdata1[2]=kap2i;
	      bdata1[3]=0;
	      bdata1[4]=0;
	      adata1[1][1]=-aa*cos(t)*k_[1]+r*sin(t)*aa*k_[2];
	      adata1[1][2]=aa*cos(t)*k_[0]-aa*sin(t)*sin(t)*k_[3];
	      adata1[1][3]=r*sin(t)*((r2+a2)*k_[3]-aa*k_[0]);
	      adata1[1][4]=a2*cos(t)*sin(t)*sin(t)*k_[1]+
		r*sin(t)*(-(r2+a2)*k_[2]);
	      adata1[2][1]=r*(-k_[1])+a2*sin(t)*cos(t)*k_[2];
	      adata1[2][2]=r*(k_[0]-aa*sin(t)*sin(t)*k_[3]);
	      adata1[2][3]=aa*sin(t)*cos(t)*((r2+a2)*k_[3]-aa*k_[0]);
	      adata1[2][4]=r*aa*sin(t)*sin(t)*k_[1]-
		aa*sin(t)*cos(t)*(-(r2+a2)*k_[2]);
	      for (j=0;j<=3;j++) {
		adata1[3][j+1]=p_[j];
		adata1[4][j+1]=po_[j];
	      }
	      ludcmp_js(adata1,4,indx,&dd);
	      lubksb_js(adata1,4,indx,bdata1);
	      
	      normf = sqrt(bdata1[2]*bdata1[2]*g_dn_ph[1][1]+
			   bdata1[3]*bdata1[3]*g_dn_ph[2][2]+
			   bdata1[4]*bdata1[4]*g_dn_ph[3][3]);
	      f_[0]=0;
	      f_[1]=bdata1[2]*sqrt(g_dn_ph[1][1])/normf;
	      f_[2]=-bdata1[3]*sqrt(g_dn_ph[2][2])/normf;
	      f_[3]=bdata1[4]*sqrt(g_dn_ph[3][3])/normf;
	      
	      psi = atan2(f_[2],f_[3]);
	      if ((iscat == 0)&&(dscat == 0)) isort = 0;
	      if ((iscat == 0)&&(dscat >= 1)) isort = 1;
	      if ((iscat == 1)&&(dscat == 0)) isort = 2;
	      if ((iscat == 1)&&(dscat >= 1)) isort = 3;
	      if ((iscat >= 2)&&(dscat == 0)) isort = 4;
	      if ((iscat >= 2)&&(dscat >= 1)) isort = 5;
	      //if (Fe_phot == 1) isort = 5;
	      //fprintf(outfile,"%d %12.4e %d %12.4e %d\n",
	      //      isort, yn[2], jth, yn[3], jph);

	      //rdsh is redshift of emitted photons as measured by 
	      //(distant) observer
	      rdsh = (v_[0]*p_[0]+v_[1]*p_[1]+v_[2]*p_[2]+v_[3]*p_[3])/pdv_em;
	      rdsh3 = rdsh*rdsh*rdsh;

	      degp = deg;
	      for (iv=0;iv<=Ne;iv++) {
		nu[iv] = nu[iv]*rdsh;
		dnu[iv] = dnu[iv]*rdsh;
		//Bnu_[iv] = Bnu_[iv]*rdsh3;
		Bnu_[iv] = Bnu_[iv];
		
		//linear scale
		if (spec_model == 1) 
		  je=(int)floor((nu[iv]-e_min)/(e_max-e_min)*Ne);
		//log scale
		if (spec_model == 2) 
		  je=(int)floor((Ne)*log(nu[iv]/e_min)/log(e_max/e_min));
		if ((je >= 0)&(je < Ne)) {
		  elo = je;
		  ehi = je+1;
		  Eph = nu[iv];
		  wlo = (nu0[ehi]-Eph)/(nu0[ehi]-nu0[elo]);
		  whi = (Eph-nu0[elo])/(nu0[ehi]-nu0[elo]);
		  if (wlo < 0) {
		    printf("%ld %ld %g %g %g\n",iv,je,nu[iv],dnu[iv],Bnu_[iv]);
		  }
		  //For line emission, the photon bundle is really spread out
		  //over a range of frequency dNU, and the dNU also gets
		  //red/blue-shifted so spectral Intensity is constant
		  
		  //Nphotons
		  //Ispec[index2(jth,je)] += 1./(-yn[4])*A_fact/dnu[elo];
		  
		  //Intensity
		  if (Fe_phot >= 0) {
		    Ispec[index2(jth,elo)] += wlo*A_fact*Bnu_[iv];
		    Ispec[index2(jth,ehi)] += whi*A_fact*Bnu_[iv];
		    Ispec_s[index3(jth,elo,isort)] += wlo*A_fact*Bnu_[iv];
		    Ispec_s[index3(jth,ehi,isort)] += whi*A_fact*Bnu_[iv];
		    Ispecr[indexrth(ir,jth,elo)] += wlo*A_fact*Bnu_[iv];
		    Ispecr[indexrth(ir,jth,ehi)] += whi*A_fact*Bnu_[iv];
		    Ispecp[indexph(jth,jph,elo)] += wlo*A_fact*Bnu_[iv];
		    Ispecp[indexph(jth,jph,ehi)] += whi*A_fact*Bnu_[iv];
		    if (isort == 1) {
		      Rspecr[indexrth(ir,jth,elo)] += wlo*A_fact*Bnu_[iv];
		      Rspecr[indexrth(ir,jth,ehi)] += whi*A_fact*Bnu_[iv];
		    }
		  }
		
		  //Iron line intensity
		  if (Fe_phot == 1) {
		    Lspec[index2(jth,elo)] += wlo*A_fact*Bnu_[iv];
		    Lspec[index2(jth,ehi)] += whi*A_fact*Bnu_[iv];
		  }

		  //Stokes Parameters
		  //Q,U same units of I_nu, so must include 1/dnu factor for
		  //spectral intensity
		  
		  //reduce deg by free-free absorption
		  //deg = degp*qnur[indexre(ir,iv)];
		  
		  Qspec[index2(jth,elo)] += wlo*A_fact*Bnu_[iv]*deg*cos(2.*psi);
		  Qspec[index2(jth,ehi)] += whi*A_fact*Bnu_[iv]*deg*cos(2.*psi);
		  Qspecr[indexrth(ir,jth,elo)] += wlo*A_fact*Bnu_[iv]*deg*cos(2.*psi);
		  Qspecr[indexrth(ir,jth,ehi)] += whi*A_fact*Bnu_[iv]*deg*cos(2.*psi);
		  Qspec_s[index3(jth,elo,isort)] += 
		    wlo*A_fact*Bnu_[iv]*deg*cos(2.*psi);
		  Qspec_s[index3(jth,ehi,isort)] += 
		    whi*A_fact*Bnu_[iv]*deg*cos(2.*psi);
		  //if (cos(yn[2]) > 0) {
		  if ((cos(yn[2]) > 0)||(TWO_SIDED >= 1)) {
		    Uspec[index2(jth,elo)] += 
		      wlo*A_fact*Bnu_[iv]*deg*sin(2.*psi);
		    Uspec[index2(jth,ehi)] += 
		      whi*A_fact*Bnu_[iv]*deg*sin(2.*psi);
		    Uspecr[indexrth(ir,jth,elo)] += 
		      wlo*A_fact*Bnu_[iv]*deg*sin(2.*psi);
		    Uspecr[indexrth(ir,jth,ehi)] += 
		      whi*A_fact*Bnu_[iv]*deg*sin(2.*psi);
		    Uspec_s[index3(jth,elo,isort)] += 
		      wlo*A_fact*Bnu_[iv]*deg*sin(2.*psi);
		    Uspec_s[index3(jth,ehi,isort)] += 
		      whi*A_fact*Bnu_[iv]*deg*sin(2.*psi);
		  }
		  //if (cos(yn[2]) < 0) {
		  if ((cos(yn[2]) < 0)&&(TWO_SIDED == 0)) {
		    Uspec[index2(jth,elo)] += 
		      wlo*A_fact*Bnu_[iv]*deg*sin(-2.*psi);
		    Uspec[index2(jth,ehi)] += 
		      whi*A_fact*Bnu_[iv]*deg*sin(-2.*psi);
		    Uspecr[indexrth(ir,jth,elo)] += 
		      wlo*A_fact*Bnu_[iv]*deg*sin(-2.*psi);
		    Uspecr[indexrth(ir,jth,ehi)] += 
		      whi*A_fact*Bnu_[iv]*deg*sin(-2.*psi);
		    Uspec_s[index3(jth,elo,isort)] += 
		      wlo*A_fact*Bnu_[iv]*deg*sin(-2.*psi);
		    Uspec_s[index3(jth,ehi,isort)] += 
		      whi*A_fact*Bnu_[iv]*deg*sin(-2.*psi);
		  }
		  if (Eph < nu_cut) 
		    Iobs[indexph(jth,jph,jt)]+= A_fact*Bnu_[iv]*dnu0[elo];
		  //if ((jth == 10)&&(jph == 0))
		  //printf("%d %d %10.4e %10.4e %10.4e %10.4e\n",
		  //   iv,elo,nu[iv],dnu0[elo],Bnu_[iv],Iobs[indexph(jth,jph,jt)]);
		}
	      
		//sort by energy bins for image files
		//linear scale
		if (spec_model == 1) 
		  je=(int)floor((nu[iv]-e_min)/(e_max-e_min)*Ne_obs);
		//log scale
		if (spec_model == 2) 
		  je=(int)floor((Ne_obs)*log(nu[iv]/e_min)/log(e_max/e_min));
		if ((je >= 0)&(je < Ne)) {
		  elo = je;
		  ehi = je+1;
		  Eph = nu[iv];
		  wlo = (nui0[ehi]-Eph)/dnui0[elo];
		  whi = (Eph-nui0[elo])/dnui0[elo];
		  //if (wlo < 0) {
		  //printf("%d %d %g %g %g\n",iv,je,nu[iv],dnu[iv],Bnu_[iv]);
		  //}
		  //NO REFLECTION CASE:
		  //if ((imageprint == 1)&&(isort == 0)) {
		  //ONLY SCATTERED PHOTONS:
		  //((imageprint == 1)&&(dscat == 0)&&(iscat >= 0)) {
		  //NORMAL CASE:
		  if ((imageprint == 1)&&(isort >= 0)) {
		    je = isort;
		    if ((ix >= 0)&&(ix <= Ni)&&(iy >= 0)&&(iy <= Ni)
			&&(fmod((int)jph,Nph_obs/4) == fmod(view_jph,Nph_obs/4))) {
		      //if (iv == 100) printf("%g %g %g %d %d %d\n",A_fact,Bnu_[iv],dnui0[elo],jph,ix,iy);
		      //if (cos(yn[2]) >= 0) {
		      if ((cos(yn[2]) >= 0)||(TWO_SIDED >= 1)) {
			spcimage[indexspci(jth,Ni-ix,Ni-iy,je)] +=
			  A_fact*Bnu_[iv]*dnui0[elo];
			spcimagex[indexspci(jth,Ni-ix,Ni-iy,je)] +=
			  A_fact*Bnu_[iv]*dnui0[elo]*deg*cos(2.*psi);
			spcimagey[indexspci(jth,Ni-ix,Ni-iy,je)] +=
			  A_fact*Bnu_[iv]*dnui0[elo]*deg*sin(2.*psi);
		      }
		      //if (cos(yn[2]) < 0) {
		      if ((cos(yn[2]) < 0)&&(TWO_SIDED == 0)) {
			spcimage[indexspci(jth,Ni-ix,iy,je)] +=
			  A_fact*Bnu_[iv]*dnui0[elo];
			spcimagex[indexspci(jth,Ni-ix,iy,je)] +=
			  A_fact*Bnu_[iv]*dnui0[elo]*deg*cos(2.*psi);
			spcimagey[indexspci(jth,Ni-ix,iy,je)] +=
			  A_fact*Bnu_[iv]*dnui0[elo]*deg*sin(-2.*psi);
		      }
		    }
		  }
		}
	      }
	      if ((jth == 1000)&&((Ni-ix)==30)&&((Ni-iy)==40)) {
		for (je=0;je<=Ne;je++) {
		  if (((je%6) == 5)||(je == Ne)) {
		    printf("%12.5e\n",nu0[je]);
		  } else {
		    printf("%12.5e ",nu0[je]);
		  }
		}
		for (je=0;je<=Ne;je++) {
		  if (((je%6) == 5)||(je == Ne)) {
		    printf("%12.5e\n",A_fact*Bnu_[je]);
		  } else {
		    printf("%12.5e ",A_fact*Bnu_[je]);
		  }
		}
	      }

	      //NO REFLECTION CASE:
	      //if ((imageprint == 1)&&(isort == 0)) {
	      //ONLY SCATTERED PHOTONS:
	      //if ((imageprint == 1)&&(dscat == 0)&&(iscat > 0)) {
	      //NORMAL CASE:
	      //printf("%g %d %d %d\n",cos(y2[2]),jph,ix,iy);
	      if ((imageprint == 1)&&(isort >= 0)) {
		if ((ix >= 0)&&(ix <= Ni)&&(iy >= 0)&&(iy <= Ni)&&(iscat >=0)
		    &&(fmod((int)jph,(Nph_obs+1)/4) == fmod(view_jph,(Nph_obs+1)/4))) {
		  //image contains energy-integrated flux, so we should
		  //not include a factor of 1/dnu
		  //printf("%g %d %d %d %d\n",cos(y2[2]),jth,jph,ix,iy);
		  //if (cos(yn[2]) >= 0) {
		  if ((cos(yn[2]) >= 0)||(TWO_SIDED >= 1)) {
		    image[indexi(jth,Ni-ix,Ni-iy)] += Iobs[indexph(jth,jph,jt)];
		    imagex[indexi(jth,Ni-ix,Ni-iy)] +=
		      Iobs[indexph(jth,jph,jt)]*deg*cos(2.*psi);
		    imagey[indexi(jth,Ni-ix,Ni-iy)] +=
		      Iobs[indexph(jth,jph,jt)]*deg*sin(2.*psi);
		  }
		  //if (cos(yn[2]) < 0) {
		  if ((cos(yn[2]) < 0)&&(TWO_SIDED == 0)) {
		    image[indexi(jth,Ni-ix,iy)] += Iobs[indexph(jth,jph,jt)];
		    imagex[indexi(jth,Ni-ix,iy)] += 
		      Iobs[indexph(jth,jph,jt)]*deg*cos(2.*psi);
		    imagey[indexi(jth,Ni-ix,iy)] += 
		      Iobs[indexph(jth,jph,jt)]*deg*(-sin(2.*psi));
		  }
		}
		if ((jth==10)&&(jph==0)) {
		  //printf("%10.4e %10.4e %10.4e %10.4e %10.4e\n",
		  // Bnu_[0],Bnu_[50],Bnu_[100],Bnu_[150],Bnu_[200]);
		  //printf("%g %g\n",A_fact,Iobs[indexph(jth,jph,jt)]);
		}
		if ((ix >= 0)&&(ix <= Ni)&&(iy >= 0)&&(iy <= Ni)
		    &&((cos(y2[2])>0)||(TWO_SIDED >= 1))&&(iscat >=0)) {
		  for (jquad=0;jquad<=3;jquad++) {
		    jphq = fmod((int)(jph+(Nph_obs+1)/4*jquad),Nph_obs+1);
		    phimage[indexphi(jth,jphq,Ni-ix,Ni-iy)] += 
		      Iobs[indexph(jth,jph,jt)];
		    phimagex[indexphi(jth,jphq,Ni-ix,Ni-iy)] +=
		      Iobs[indexph(jth,jph,jt)]*deg*cos(2.*psi);
		    phimagey[indexphi(jth,jphq,Ni-ix,Ni-iy)] +=
		      Iobs[indexph(jth,jph,jt)]*deg*sin(2.*psi);
		  }
		}
	      }
	      Nescape++;
	      tau_tt = tau_tt+tau_tot;
	    }
	    //printf("%d %12.5g %12.5g %12.5g %12.5g %12.5g\n",
	    // iscat,tau_tot,l_tot/(1.45e6),yn[1],yn[4],yn[4]-y0[4]);
	  }
	  //printf("%d %d %d %d %g %g %g\n",
	  // it,ip,(int)Nescape,steps,A_fact,yn[1],yn[4]);
	  Nscat+=iscat;
	  Nscat+=dscat;
	} //ip=0..N
      } //it=0..N
      if ((myid == 0)&&(iph == 0)&&(fmod(ir,5)==0)&&(ibottom>=0)) 
	printf("%d %d %g %g %g %g %g %g\n",
	       (int)RUN_ID,ibottom,rr[ir],Nescape,Ndisk,Nscat,
	       emtop_elf_ik[indexelf(ir,iph,3,1)],tau_tt/((N+1)*(N+1)));
      /*
      if ((myid == 0)&&(iph == 0)&&(fmod(ir,10)==0)&&(ibottom==1)) 
	printf("%d %d %g %g %g %g %g %g\n",
	       (int)RUN_ID,iph,rr[ir],Nescape,Ndisk,
	       embot_elf_ik[indexelf(ir,iph,3,1)],
	       cos(embot_ik[indexr(ir,iph)]),Nscat/((N+1)*(N+1)));
      */
    } //iph=0..Nph
    //pass_data(myid,numprocs,tag,comm,Ispec,Lspec,Rspec,I_r);
  } //ir=0..Nr
  } //ibottom 
  //schem_ray2
  //fclose(outfile);

  if (myid == 0) { 
    file_ik = 0;
    file_id = RUN_ID;
    //    if (OPT_THIN == 1) file_id++;
    file_ik = (file_id-fmod(file_id,1000))/1000;
    file_id = file_id-file_ik*1000;
    file_ih = (file_id-fmod(file_id,100))/100;
    file_id = file_id-file_ih*100;
    file_iu = fmod(file_id,10);
    file_id = (file_id-fmod(file_id,10))/10;
    strcpy(fname1,"data/scat_spec.0000.dat");
    fname1[15]=48+file_ik;
    fname1[16]=48+file_ih;
    fname1[17]=48+file_id;
    fname1[18]=48+file_iu;
    strcpy(fname2,"data/scat_line.0000.dat");
    fname2[15]=48+file_ik;
    fname2[16]=48+file_ih;
    fname2[17]=48+file_id;
    fname2[18]=48+file_iu;
    strcpy(fname3,"data/scat_disk.0000.dat");
    fname3[15]=48+file_ik;
    fname3[16]=48+file_ih;
    fname3[17]=48+file_id;
    fname3[18]=48+file_iu;
    strcpy(fname4,"data/scat_imag.0000.dat");
    fname4[15]=48+file_ik;
    fname4[16]=48+file_ih;
    fname4[17]=48+file_id;
    fname4[18]=48+file_iu;
    strcpy(fname5,"data/scat_ipol.0000.dat");
    fname5[15]=48+file_ik;
    fname5[16]=48+file_ih;
    fname5[17]=48+file_id;
    fname5[18]=48+file_iu;
    strcpy(fname6,"data/scat_spcr.0000.dat");
    fname6[15]=48+file_ik;
    fname6[16]=48+file_ih;
    fname6[17]=48+file_id;
    fname6[18]=48+file_iu;
    strcpy(fname7,"data/scat_ithp.0000.dat");
    fname7[15]=48+file_ik;
    fname7[16]=48+file_ih;
    fname7[17]=48+file_id;
    fname7[18]=48+file_iu;
    strcpy(fname8,"data/scat_cpow.0000.dat");
    fname8[15]=48+file_ik;
    fname8[16]=48+file_ih;
    fname8[17]=48+file_id;
    fname8[18]=48+file_iu;
    strcpy(fname9,"data/scat_spcp.0000.dat");
    fname9[15]=48+file_ik;
    fname9[16]=48+file_ih;
    fname9[17]=48+file_id;
    fname9[18]=48+file_iu;

    outfile = fopen(fname1, "w");
    outfile2 = fopen(fname2, "w");
    for (jth=0;jth<=Nth_obs;jth++) {
      for (je=0;je<=Ne;je++) {
	if (((je%6) == 5)||(je == Ne)) {
	  fprintf(outfile,"%12.5e\n",Ispec[index2(jth,je)]);
	  fprintf(outfile2,"%12.5e\n",Lspec[index2(jth,je)]);
	} else {
	  fprintf(outfile,"%12.5e ",Ispec[index2(jth,je)]);
	  fprintf(outfile2,"%12.5e ",Lspec[index2(jth,je)]);
	}
      }
    }
    fclose(outfile2);
    fprintf(outfile,"\n");
    for (jth=0;jth<=Nth_obs;jth++) {
      for (je=0;je<=Ne;je++) {
	if (((je%6) == 5)||(je == Ne)) {
	  fprintf(outfile,"%12.5e\n",Qspec[index2(jth,je)]);
	} else {
	  fprintf(outfile,"%12.5e ",Qspec[index2(jth,je)]);
	}
      }
    }
    fprintf(outfile,"\n");
    for (jth=0;jth<=Nth_obs;jth++) {
      for (je=0;je<=Ne;je++) {
	if (((je%6) == 5)||(je == Ne)) {
	  fprintf(outfile,"%12.5e\n",Uspec[index2(jth,je)]);
	} else {
	  fprintf(outfile,"%12.5e ",Uspec[index2(jth,je)]);
	}
      }
    }
    
    for (isort=0;isort<=5;isort++) {
      fprintf(outfile,"\n");
      for (jth=0;jth<=Nth_obs;jth++) {
	for (je=0;je<=Ne;je++) {
	  if (((je%6) == 5)||(je == Ne)) {
	    fprintf(outfile,"%12.5e\n",Ispec_s[index3(jth,je,isort)]);
	  } else {
	    fprintf(outfile,"%12.5e ",Ispec_s[index3(jth,je,isort)]);
	  }
	}
      }
      fprintf(outfile,"\n");
      for (jth=0;jth<=Nth_obs;jth++) {
	for (je=0;je<=Ne;je++) {
	  if (((je%6) == 5)||(je == Ne)) {
	    fprintf(outfile,"%12.5e\n",Qspec_s[index3(jth,je,isort)]);
	  } else {
	    fprintf(outfile,"%12.5e ",Qspec_s[index3(jth,je,isort)]);
	  }
	}
      }
      fprintf(outfile,"\n");
      for (jth=0;jth<=Nth_obs;jth++) {
	for (je=0;je<=Ne;je++) {
	  if (((je%6) == 5)||(je == Ne)) {
	    fprintf(outfile,"%12.5e\n",Uspec_s[index3(jth,je,isort)]);
	  } else {
	    fprintf(outfile,"%12.5e ",Uspec_s[index3(jth,je,isort)]);
	  }
	}
      }
    }
    fclose(outfile);
    outfile = fopen(fname3, "w");
    outfile2 = fopen(fname6, "w");
    for (jr=0;jr<=Nr;jr++) {
      for (je=0;je<=Ne;je++) {
	//printf("%d %d %g\n",jr,je,Rspec[indexre(jr,je)]);

	if (((je%6) == 5)||(je == Ne)) {
	  fprintf(outfile,"%12.5e\n",Rspec[indexre(jr,je)]);
	} else {
	  fprintf(outfile,"%12.5e ",Rspec[indexre(jr,je)]);
	}
      }

      for (jth=0;jth<=Nth_obs;jth++) {
	for (je=0;je<=Ne;je++) {
	  if (((je%6) == 5)||(je == Ne)) {
	    fprintf(outfile2,"%12.5e\n",Ispecr[indexrth(jr,jth,je)]);
	  } else {
	    fprintf(outfile2,"%12.5e ",Ispecr[indexrth(jr,jth,je)]);
	  }
	}
      }
      fprintf(outfile2,"\n");
    }
    fclose(outfile);

    for (jr=0;jr<=Nr;jr++) {
      fprintf(outfile2,"%12.5e %12.5e %12.5e\n",
	      rr[jr],L_factr[jr],G_factr[jr]);
    }
    fprintf(outfile2,"\n");

    for (jr=0;jr<=Nr;jr++) {
      for (jth=0;jth<=Nth_obs;jth++) {
	for (je=0;je<=Ne;je++) {
	  if (((je%6) == 5)||(je == Ne)) {
	    fprintf(outfile2,"%12.5e\n",Qspecr[indexrth(jr,jth,je)]);
	  } else {
	    fprintf(outfile2,"%12.5e ",Qspecr[indexrth(jr,jth,je)]);
	  }
	}
      }
      fprintf(outfile2,"\n");
    }
    fprintf(outfile2,"\n");

    for (jr=0;jr<=Nr;jr++) {
      for (jth=0;jth<=Nth_obs;jth++) {
	for (je=0;je<=Ne;je++) {
	  if (((je%6) == 5)||(je == Ne)) {
	    fprintf(outfile2,"%12.5e\n",Uspecr[indexrth(jr,jth,je)]);
	  } else {
	    fprintf(outfile2,"%12.5e ",Uspecr[indexrth(jr,jth,je)]);
	  }
	}
      }
      fprintf(outfile2,"\n");
    }
    fprintf(outfile2,"\n");

    for (jr=0;jr<=Nr;jr++) {
      for (jth=0;jth<=Nth_obs;jth++) {
	for (je=0;je<=Ne;je++) {
	  if (((je%6) == 5)||(je == Ne)) {
	    fprintf(outfile2,"%12.5e\n",Rspecr[indexrth(jr,jth,je)]);
	  } else {
	    fprintf(outfile2,"%12.5e ",Rspecr[indexrth(jr,jth,je)]);
	  }
	}
      }
      fprintf(outfile2,"\n");
    }
    for (jr=0;jr<=Nr;jr++) {
      for (je=0;je<=Ne;je++) {
	if (((je%6) == 5)||(je == Ne)) {
	  fprintf(outfile2,"%12.5e\n",Cspecr[indexre(jr,je)]);
	} else {
	  fprintf(outfile2,"%12.5e ",Cspecr[indexre(jr,je)]);
	}
      }
    }
    fclose(outfile2);

    outfile = fopen(fname4, "w");
    for (it=0;it<=Nth_obs;it++) {
      for (ix=0;ix<=Ni;ix++) {
	for (iy=0;iy<=Ni;iy++) {
	  if (((iy%6) == 5)||(iy == Ni)) {
	    fprintf(outfile,"%12.5e\n",image[indexi(it,ix,iy)]);
	  } else {
	    fprintf(outfile,"%12.5e ",image[indexi(it,ix,iy)]);
	  }
	}
      }
      fprintf(outfile,"\n");
    }

    for (it=0;it<=Nth_obs;it++) {
      for (ix=0;ix<=Ni;ix++) {
	for (iy=0;iy<=Ni;iy++) {
	  for (je=0;je<=Ne_obs;je++) {
	    if (((je%6) == 5)||(je == Ne_obs)) {
	      fprintf(outfile,"%12.5e\n",spcimage[indexspci(it,ix,iy,je)]);
	    } else {
	      fprintf(outfile,"%12.5e ",spcimage[indexspci(it,ix,iy,je)]);
	    }
	  }
	}
      }
      fprintf(outfile,"\n");
    }
    fclose(outfile);

    outfile = fopen(fname5, "w");
    for (it=0;it<=Nth_obs;it++) {
      for (ix=0;ix<=Ni;ix++) {
	for (iy=0;iy<=Ni;iy++) {
	  if (((iy%6) == 5)||(iy == Ni)) {
	    fprintf(outfile,"%12.5e\n",imagex[indexi(it,ix,iy)]);
	  } else {
	    fprintf(outfile,"%12.5e ",imagex[indexi(it,ix,iy)]);
	  }
	}
      }
      fprintf(outfile,"\n");
    }
    for (it=0;it<=Nth_obs;it++) {
      for (ix=0;ix<=Ni;ix++) {
	for (iy=0;iy<=Ni;iy++) {
	  if (((iy%6) == 5)||(iy == Ni)) {
	    fprintf(outfile,"%12.5e\n",imagey[indexi(it,ix,iy)]);
	  } else {
	    fprintf(outfile,"%12.5e ",imagey[indexi(it,ix,iy)]);
	  }
	}
      }
      fprintf(outfile,"\n");
    }

    for (it=0;it<=Nth_obs;it++) {
      for (ix=0;ix<=Ni;ix++) {
	for (iy=0;iy<=Ni;iy++) {
	  for (je=0;je<=Ne_obs;je++) {
	    if (((je%6) == 5)||(je == Ne_obs)) {
	      fprintf(outfile,"%12.5e\n",spcimagex[indexspci(it,ix,iy,je)]);
	    } else {
	      fprintf(outfile,"%12.5e ",spcimagex[indexspci(it,ix,iy,je)]);
	    }
	  }
	}
      }
      fprintf(outfile,"\n");
    }
    for (it=0;it<=Nth_obs;it++) {
      for (ix=0;ix<=Ni;ix++) {
	for (iy=0;iy<=Ni;iy++) {
	  for (je=0;je<=Ne_obs;je++) {
	    if (((je%6) == 5)||(je == Ne_obs)) {
	      fprintf(outfile,"%12.5e\n",spcimagey[indexspci(it,ix,iy,je)]);
	    } else {
	      fprintf(outfile,"%12.5e ",spcimagey[indexspci(it,ix,iy,je)]);
	    }
	  }
	}
      }
      fprintf(outfile,"\n");
    }
    fclose(outfile);

    outfile = fopen(fname8, "w");
    for (ir=0;ir<=Nr;ir++) {
      for (ith=0;ith<=Nth;ith++) {
	for (iph=0;iph<=Nph;iph++) {
	  if (((iph%6) == 5)||(iph == Nph)) {
	    fprintf(outfile,"%12.5e\n",corpow_ijk[indexijk(ir,ith,iph)]);
	  } else {
	    fprintf(outfile,"%12.5e ",corpow_ijk[indexijk(ir,ith,iph)]);
	  }
	}
      }
      fprintf(outfile,"\n");
    }
    fclose(outfile);

    outfile = fopen(fname9, "w");
    for (jth=0;jth<=Nth_obs;jth++) {
      for (jph=0;jph<=Nph_obs;jph++) {
	for (je=0;je<=Ne;je++) {
	  if (((je%6) == 5)||(je == Ne)) {
	    fprintf(outfile,"%12.5e\n",Ispecp[indexph(jth,jph,je)]);
	  } else {
	    fprintf(outfile,"%12.5e ",Ispecp[indexph(jth,jph,je)]);
	  }
	}
      }
      fprintf(outfile,"\n");
    }
    fclose(outfile);

    /*
    for (it=0;it<=Nth_obs;it++) {
      for (ix=0;ix<=Ni;ix++) {
	for (iy=0;iy<=Ni;iy++) {
	  if (((iy%6) == 5)||(iy == Ni)) {
	    fprintf(outfile,"%12.5e\n",imagey[indexi(it,ix,iy)]);
	  } else {
	    fprintf(outfile,"%12.5e ",imagey[indexi(it,ix,iy)]);
	  }
	}
      }
      fprintf(outfile,"\n");
    }
    */

    outfile = fopen(fname7, "w");
    for (it=0;it<=Nth_obs;it++) {
      for (ip=0;ip<=Nph_obs;ip++) {
	for (ix=0;ix<=Ni;ix++) {
	  for (iy=0;iy<=Ni;iy++) {
	    if (((iy%6) == 5)||(iy == Ni)) {
	      fprintf(outfile,"%12.5e\n",phimage[indexphi(it,ip,ix,iy)]);
	    } else {
	      fprintf(outfile,"%12.5e ",phimage[indexphi(it,ip,ix,iy)]);
	    }
	  }
	}
	fprintf(outfile,"\n");
      }
      fprintf(outfile,"\n");
    }
    for (it=0;it<=Nth_obs;it++) {
      for (ip=0;ip<=Nph_obs;ip++) {
	for (ix=0;ix<=Ni;ix++) {
	  for (iy=0;iy<=Ni;iy++) {
	    if (((iy%6) == 5)||(iy == Ni)) {
	      fprintf(outfile,"%12.5e\n",phimagex[indexphi(it,ip,ix,iy)]);
	    } else {
	      fprintf(outfile,"%12.5e ",phimagex[indexphi(it,ip,ix,iy)]);
	    }
	  }
	}
	fprintf(outfile,"\n");
      }
      fprintf(outfile,"\n");
    }
    for (it=0;it<=Nth_obs;it++) {
      for (ip=0;ip<=Nph_obs;ip++) {
	for (ix=0;ix<=Ni;ix++) {
	  for (iy=0;iy<=Ni;iy++) {
	    if (((iy%6) == 5)||(iy == Ni)) {
	      fprintf(outfile,"%12.5e\n",phimagey[indexphi(it,ip,ix,iy)]);
	    } else {
	      fprintf(outfile,"%12.5e ",phimagey[indexphi(it,ip,ix,iy)]);
	    }
	  }
	}
	fprintf(outfile,"\n");
      }
      fprintf(outfile,"\n");
    }
    fclose(outfile);
  }
  wrt_time();
  //MPI_Finalize();
  return(0);
}

