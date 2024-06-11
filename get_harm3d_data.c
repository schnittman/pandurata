#include "panhead.h"

void get_harm3d_data(double rr[], double tt[], double pp[],
		     double rho_ijk[], double T_ijk[], double bb_ijk[], 
		     double tau_ijk[], double ut_ijk[], double ur_ijk[], 
		     double uz_ijk[], double up_ijk[],
		     int diskbody_ik[], double sigtau_ik[], double Tdisk_ik[],
		     double emtop_ik[], double embot_ik[], 
		     double reftop_ik[], double refbot_ik[])
{
  double blnk1, blnk2, blnk3;
  float *rho_flt, *T_flt, *bb_flt, *tau_flt, *ut_flt, *ur_flt, *uz_flt, 
    *up_flt, *photo_flt;
  long i,j,k,ik,ih,id,iu;
  FILE *grid_file, *rho_file, *T_file, *bb_file, *tau_file, *ut_file, 
    *ur_file, *uz_file, *up_file, *photo_file;
  char fname[100];

  //build filename string
  ik = 0;
  id = (long)RUN_ID;
  ik = (id-fmod(id,1000))/1000;
  id = id-ik*1000;
  ih = (id-fmod(id,100))/100;
  id = id-ih*100;
  iu = fmod(id,10);
  id = (id-fmod(id,10))/10;
  printf("check\n");

  strcpy(fname,"data/gr_0000.dat");
  fname[8]=48+ik;
  fname[9]=48+ih;
  fname[10]=48+id;
  fname[11]=48+iu;
  printf("%s\n",fname);
  grid_file = fopen(fname, "r");
  fscanf(grid_file,"%lf %lf %lf\n", &blnk1, &blnk2, &blnk3);
  for (i=0;i<=Nr;i++) {
    if (((i%6) == 5)||(i == Nr)) {
      fscanf(grid_file,"%lf\n", &blnk1);
    } else {
      fscanf(grid_file,"%lf", &blnk1);
    }
    rr[i]=blnk1;
  }
  for (i=0;i<=Nth;i++) {
    if (((i%6) == 5)||(i == Nth)) {
      fscanf(grid_file,"%lf\n", &blnk1);
    } else {
      fscanf(grid_file,"%lf", &blnk1);
    }
    tt[i]=blnk1;
  }
  for (i=0;i<=Nph;i++) {
    if (((i%6) == 5)||(i == Nph)) {
      fscanf(grid_file,"%lf\n", &blnk1);
    } else {
      fscanf(grid_file,"%lf", &blnk1);
    }
    pp[i]=blnk1;
  }
  fclose(grid_file);

  strcpy(fname,"data/ph_0000.dat");
  fname[8]=48+ik;
  fname[9]=48+ih;
  fname[10]=48+id;
  fname[11]=48+iu;
  photo_file = fopen(fname, "r");
  //read in diskbody: index of optical depth in (r,theta)
  for (k=0;k<=Nph;k++) {
    for (i=0;i<=Nr;i++) {
      if (((i%6) == 5)||(i == Nr)) {
	fscanf(photo_file,"%lf\n", &blnk1);
      } else {
	fscanf(photo_file,"%lf", &blnk1);
      }
      diskbody_ik[indexr(i,k)]=(int)blnk1;
    }
  }      
  //read in sigtau_es
  for (k=0;k<=Nph;k++) {
    for (i=0;i<=Nr;i++) {
      if (((i%6) == 5)||(i == Nr)) {
	fscanf(photo_file,"%lf\n", &blnk1);
      } else {
	fscanf(photo_file,"%lf", &blnk1);
      }
      sigtau_ik[indexr(i,k)]=blnk1;
    }
  }      
  //read in Tdisk_ik: surface temperature of disk
  for (k=0;k<=Nph;k++) {
    for (i=0;i<=Nr;i++) {
      if (((i%6) == 5)||(i == Nr)) {
	fscanf(photo_file,"%lf\n", &blnk1);
      } else {
	fscanf(photo_file,"%lf", &blnk1);
      }
      Tdisk_ik[indexr(i,k)]=blnk1;
      //scale up temp to match spectra to 0.1 Ledd
      //Tdisk_ik[indexr(i,k)]=Tdisk_ik[indexr(i,k)]*pow(0.1/L_Edd,0.25);
    }
  }      
  //read in emtop_ik: theta coordinates of emission surface
  for (k=0;k<=Nph;k++) {
    for (i=0;i<=Nr;i++) {
      if (((i%6) == 5)||(i == Nr)) {
	fscanf(photo_file,"%lf\n", &blnk1);
      } else {
	fscanf(photo_file,"%lf", &blnk1);
      }
      emtop_ik[indexr(i,k)]=blnk1;
    }
  }      
  //read in embot_ik: theta coordinates of emission surface
  for (k=0;k<=Nph;k++) {
    for (i=0;i<=Nr;i++) {
      if (((i%6) == 5)||(i == Nr)) {
	fscanf(photo_file,"%lf\n", &blnk1);
      } else {
	fscanf(photo_file,"%lf", &blnk1);
      }
      embot_ik[indexr(i,k)]=blnk1;
    }
  }      
  //read in reftop_ik: theta coordinates of reflection surface
  for (k=0;k<=Nph;k++) {
    for (i=0;i<=Nr;i++) {
      if (((i%6) == 5)||(i == Nr)) {
	fscanf(photo_file,"%lf\n", &blnk1);
      } else {
	fscanf(photo_file,"%lf", &blnk1);
      }
      reftop_ik[indexr(i,k)]=blnk1;
    }
  }      
  //read in refbot_ik: theta coordinates of reflection surface
  for (k=0;k<=Nph;k++) {
    for (i=0;i<=Nr;i++) {
      if (((i%6) == 5)||(i == Nr)) {
	fscanf(photo_file,"%lf\n", &blnk1);
      } else {
	fscanf(photo_file,"%lf", &blnk1);
      }
      refbot_ik[indexr(i,k)]=blnk1;
    }
  }      
  fclose(photo_file);

  rho_flt = (float *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(float));
  T_flt = (float *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(float));
  bb_flt = (float *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(float));
  tau_flt = (float *)malloc((Nr+1)*(Nth+2)*(Nph+1)*sizeof(float));
  ut_flt = (float *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(float));
  ur_flt = (float *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(float));
  uz_flt = (float *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(float));
  up_flt = (float *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(float));

  strcpy(fname,"data/rh_0000.bin");
  fname[8]=48+ik;
  fname[9]=48+ih;
  fname[10]=48+id;
  fname[11]=48+iu;
  rho_file = fopen(fname, "rb");
  fread(rho_flt,sizeof(float),(Nr+1)*(Nth+1)*(Nph+1),rho_file);
  fclose(rho_file);

  strcpy(fname,"data/te_0000.bin");
  fname[8]=48+ik;
  fname[9]=48+ih;
  fname[10]=48+id;
  fname[11]=48+iu;
  T_file = fopen(fname, "rb");
  fread(T_flt,sizeof(float),(Nr+1)*(Nth+1)*(Nph+1),T_file);
  fclose(T_file);

  /*
  strcpy(fname,"data/bb_0000.bin");
  fname[8]=48+ik;
  fname[9]=48+ih;
  fname[10]=48+id;
  fname[11]=48+iu;
  bb_file = fopen(fname, "rb");
  fread(bb_flt,sizeof(float),(Nr+1)*(Nth+1)*(Nph+1),T_file);
  fclose(bb_file);
  */

  strcpy(fname,"data/ta_0000.bin");
  fname[8]=48+ik;
  fname[9]=48+ih;
  fname[10]=48+id;
  fname[11]=48+iu;
  tau_file = fopen(fname, "rb");
  fread(tau_flt,sizeof(float),(Nr+1)*(Nth+2)*(Nph+1),tau_file);
  fclose(tau_file);

  strcpy(fname,"data/u0_0000.bin");
  fname[8]=48+ik;
  fname[9]=48+ih;
  fname[10]=48+id;
  fname[11]=48+iu;
  ut_file = fopen(fname, "rb");
  fread(ut_flt,sizeof(float),(Nr+1)*(Nth+1)*(Nph+1),ut_file);
  fclose(ut_file);

  strcpy(fname,"data/u1_0000.bin");
  fname[8]=48+ik;
  fname[9]=48+ih;
  fname[10]=48+id;
  fname[11]=48+iu;
  ur_file = fopen(fname, "rb");
  fread(ur_flt,sizeof(float),(Nr+1)*(Nth+1)*(Nph+1),ur_file);
  fclose(ur_file);

  strcpy(fname,"data/u2_0000.bin");
  fname[8]=48+ik;
  fname[9]=48+ih;
  fname[10]=48+id;
  fname[11]=48+iu;
  uz_file = fopen(fname, "rb");
  fread(uz_flt,sizeof(float),(Nr+1)*(Nth+1)*(Nph+1),uz_file);
  fclose(uz_file);

  strcpy(fname,"data/u3_0000.bin");
  fname[8]=48+ik;
  fname[9]=48+ih;
  fname[10]=48+id;
  fname[11]=48+iu;
  up_file = fopen(fname, "rb");
  fread(up_flt,sizeof(float),(Nr+1)*(Nth+1)*(Nph+1),up_file);
  fclose(up_file);

  for (i=0;i<=Nr;i++) {
    for (j=0;j<=Nth;j++) {
      for (k=0;k<=Nph;k++) {
	rho_ijk[indexijk(i,j,k)]=(double)rho_flt[indexijk(i,j,k)];
	T_ijk[indexijk(i,j,k)]=(double)T_flt[indexijk(i,j,k)];
	bb_ijk[indexijk(i,j,k)]=(double)bb_flt[indexijk(i,j,k)];
	ut_ijk[indexijk(i,j,k)]=(double)ut_flt[indexijk(i,j,k)];
	ur_ijk[indexijk(i,j,k)]=(double)ur_flt[indexijk(i,j,k)];
	uz_ijk[indexijk(i,j,k)]=(double)uz_flt[indexijk(i,j,k)];
	up_ijk[indexijk(i,j,k)]=(double)up_flt[indexijk(i,j,k)];
      }
    }
  }
  for (i=0;i<=Nr;i++) {
    for (j=0;j<=Nth+1;j++) {
      for (k=0;k<=Nph;k++) {
	tau_ijk[indexthijk(i,j,k)]=(double)tau_flt[indexthijk(i,j,k)];
      }
    }
  }
  return;
}

