#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define Nr 1823
//#define Nr 191
//#define Nr 203
#define Nr 179 //ThinHR a=0.0
//#define Nr 191 //ThinHR a=0.5
//#define Nr 203 //ThinHR a=0.9
//#define Nr 209 //ThinHR a=0.99
//#define Nr 359 //ThinHR+NT a=0.0
//#define Nr 383 //ThinHR+NT a=0.5
//#define Nr 407 //ThinHR+NT a=0.9
//#define Nr 419 //ThinHR+NT a=0.99
#define Nth 159
//#define Nth 191
#define Nph 63
#define Nstart 1250
#define Nstop 1250
#define Nstep 2

#define index2(a,b) ((Nfreq+1)*(b)+(a))
#define index3(a,b,c) (8*(N+1)*(c)+8*(b)+(a))
#define index4(a,b,c,d) (4*(Nlat+1)*(N+1)*(d)+4*(Nlat+1)*(c)+4*(b)+(a))
#define indexdex(b,c,d) ((Nlat+1)*(N+1)*(d)+(Nlat+1)*(c)+(b))
#define indexijk(a,b,c) ((Nph+1)*(Nth+1)*(a)+(Nph+1)*(b)+(c))
#define indexthijk(a,b,c) ((Nph+1)*(Nth+2)*(a)+(Nph+1)*(b)+(c))
#define movdex(a,b) ((N+1)*(b)+(a))

int main(void)
{
  float *rho_ijk,*T_ijk,*bb_ijk,*tau_ijk,*Ut_ijk,*Ur_ijk,*Uz_ijk,*Up_ijk,
    *phtop_ik,*phbot_ik,z1,z2,z3,z4,z5,z6;
  long it,i_start,i_stop,i,j,k,
    ilat,ik,id,ih,iu,ilo,ihi,jlo,jhi,klo,khi,i10,j10,k10;
  FILE *infile,*outfile,
    *rho_file,*T_file,*bb_file,*Ut_file,*Ur_file,*Uz_file,*Up_file,*tau_file;
  char fname[14];

  T_ijk = (float *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(float));
  rho_ijk = (float *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(float));
  Ut_ijk = (float *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(float));
  Ur_ijk = (float *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(float));
  Uz_ijk = (float *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(float));
  Up_ijk = (float *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(float));
  tau_ijk = (float *)malloc((Nr+1)*(Nth+2)*(Nph+1)*sizeof(float));
  bb_ijk = (float *)malloc((Nr+1)*(Nth+1)*(Nph+1)*sizeof(float));

  for (it=Nstart;it<=Nstop;it+=Nstep) {
    //build filename string
    ik = 0;
    id = it;
    ik = (id-fmod(id,1000))/1000;
    id = id-ik*1000;
    ih = (id-fmod(id,100))/100;
    id = id-ih*100;
    iu = fmod(id,10);
    id = (id-fmod(id,10))/10;

    printf("%d %d %d %d %d\n",it,ik,ih,id,iu);

    strcpy(fname,"rh_0000.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    rho_file = fopen(fname, "r");
    strcpy(fname,"te_0000.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    T_file = fopen(fname, "r");
    strcpy(fname,"u0_0000.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    Ut_file = fopen(fname, "r");
    strcpy(fname,"u1_0000.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    Ur_file = fopen(fname, "r");
    strcpy(fname,"u2_0000.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    Uz_file = fopen(fname, "r");
    strcpy(fname,"u3_0000.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    Up_file = fopen(fname, "r");
    
    /*
    strcpy(fname,"bb_0000.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    bb_file = fopen(fname, "r");
    */
    /*
    rho_file = fopen("rho_data_new.dat","r");
    T_file = fopen("T_data_new.dat","r");
    Ut_file = fopen("u0_data_new.dat","r");
    Ur_file = fopen("u1_data_new.dat","r");
    Uz_file = fopen("u2_data_new.dat","r");
    Up_file = fopen("u3_data_new.dat","r");
    tau_file = fopen("tau_data_new.dat","r");
    rho_file = fopen("rho_data_fake.dat","r");
    T_file = fopen("T_data_fake.dat","r");
    Ut_file = fopen("u0_data_fake.dat","r");
    Ur_file = fopen("u1_data_fake.dat","r");
    Uz_file = fopen("u2_data_fake.dat","r");
    Up_file = fopen("u3_data_fake.dat","r");
    tau_file = fopen("tau_data_fake.dat","r");
    */
    for (k=0;k<=Nph;k++) {
      for (j=0;j<=Nth;j++) {
	for (i=0;i<=Nr;i+=6) {
	  fscanf(rho_file,"%f %f %f %f %f %f\n",&z1,&z2,&z3,&z4,&z5,&z6);
	  rho_ijk[indexijk(i,j,k)]=z1;
	  rho_ijk[indexijk(i+1,j,k)]=z2;
	  rho_ijk[indexijk(i+2,j,k)]=z3;
	  rho_ijk[indexijk(i+3,j,k)]=z4;
	  rho_ijk[indexijk(i+4,j,k)]=z5;
	  rho_ijk[indexijk(i+5,j,k)]=z6;
	  fscanf(T_file,"%f %f %f %f %f %f\n",&z1,&z2,&z3,&z4,&z5,&z6);
	  T_ijk[indexijk(i,j,k)]=z1;
	  T_ijk[indexijk(i+1,j,k)]=z2;
	  T_ijk[indexijk(i+2,j,k)]=z3;
	  T_ijk[indexijk(i+3,j,k)]=z4;
	  T_ijk[indexijk(i+4,j,k)]=z5;
	  T_ijk[indexijk(i+5,j,k)]=z6;
	  fscanf(Ut_file,"%f %f %f %f %f %f\n",&z1,&z2,&z3,&z4,&z5,&z6);
	  Ut_ijk[indexijk(i,j,k)]=z1;
	  Ut_ijk[indexijk(i+1,j,k)]=z2;
	  Ut_ijk[indexijk(i+2,j,k)]=z3;
	  Ut_ijk[indexijk(i+3,j,k)]=z4;
	  Ut_ijk[indexijk(i+4,j,k)]=z5;
	  Ut_ijk[indexijk(i+5,j,k)]=z6;
	  fscanf(Ur_file,"%f %f %f %f %f %f\n",&z1,&z2,&z3,&z4,&z5,&z6);
	  Ur_ijk[indexijk(i,j,k)]=z1;
	  Ur_ijk[indexijk(i+1,j,k)]=z2;
	  Ur_ijk[indexijk(i+2,j,k)]=z3;
	  Ur_ijk[indexijk(i+3,j,k)]=z4;
	  Ur_ijk[indexijk(i+4,j,k)]=z5;
	  Ur_ijk[indexijk(i+5,j,k)]=z6;
	  fscanf(Uz_file,"%f %f %f %f %f %f\n",&z1,&z2,&z3,&z4,&z5,&z6);
	  Uz_ijk[indexijk(i,j,k)]=z1;
	  Uz_ijk[indexijk(i+1,j,k)]=z2;
	  Uz_ijk[indexijk(i+2,j,k)]=z3;
	  Uz_ijk[indexijk(i+3,j,k)]=z4;
	  Uz_ijk[indexijk(i+4,j,k)]=z5;
	  Uz_ijk[indexijk(i+5,j,k)]=z6;
	  fscanf(Up_file,"%f %f %f %f %f %f\n",&z1,&z2,&z3,&z4,&z5,&z6);
	  Up_ijk[indexijk(i,j,k)]=z1;
	  Up_ijk[indexijk(i+1,j,k)]=z2;
	  Up_ijk[indexijk(i+2,j,k)]=z3;
	  Up_ijk[indexijk(i+3,j,k)]=z4;
	  Up_ijk[indexijk(i+4,j,k)]=z5;
	  Up_ijk[indexijk(i+5,j,k)]=z6;
	  /*
	  fscanf(bb_file,"%f %f %f %f %f %f\n",&z1,&z2,&z3,&z4,&z5,&z6);
	  bb_ijk[indexijk(i,j,k)]=z1;
	  bb_ijk[indexijk(i+1,j,k)]=z2;
	  bb_ijk[indexijk(i+2,j,k)]=z3;
	  bb_ijk[indexijk(i+3,j,k)]=z4;
	  bb_ijk[indexijk(i+4,j,k)]=z5;
	  bb_ijk[indexijk(i+5,j,k)]=z6;
	  */
	}
      }
    } 
    fclose(rho_file);
    fclose(T_file);
    fclose(Ut_file);
    fclose(Ur_file);
    fclose(Uz_file);
    fclose(Up_file);
    //fclose(bb_file);
    strcpy(fname,"ta_0000.dat");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    tau_file = fopen(fname, "r");
    for (k=0;k<=Nph;k++) {
      for (j=0;j<=Nth+1;j++) {
	for (i=0;i<=Nr;i+=6) {
	  fscanf(tau_file,"%f %f %f %f %f %f\n",&z1,&z2,&z3,&z4,&z5,&z6);
	  tau_ijk[indexthijk(i,j,k)]=z1;
	  tau_ijk[indexthijk(i+1,j,k)]=z2;
	  tau_ijk[indexthijk(i+2,j,k)]=z3;
	  tau_ijk[indexthijk(i+3,j,k)]=z4;
	  tau_ijk[indexthijk(i+4,j,k)]=z5;
	  tau_ijk[indexthijk(i+5,j,k)]=z6;
	}
      }
    } 
    fclose(tau_file);

    strcpy(fname,"rh_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(rho_ijk,sizeof(float),(Nr+1)*(Nth+1)*(Nph+1),outfile);
    fclose(outfile);

    strcpy(fname,"te_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(T_ijk,sizeof(float),(Nr+1)*(Nth+1)*(Nph+1),outfile);
    fclose(outfile);
    
    strcpy(fname,"u0_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(Ut_ijk,sizeof(float),(Nr+1)*(Nth+1)*(Nph+1),outfile);
    fclose(outfile);

    strcpy(fname,"u1_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(Ur_ijk,sizeof(float),(Nr+1)*(Nth+1)*(Nph+1),outfile);
    fclose(outfile);

    strcpy(fname,"u2_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(Uz_ijk,sizeof(float),(Nr+1)*(Nth+1)*(Nph+1),outfile);
    fclose(outfile);

    strcpy(fname,"u3_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(Up_ijk,sizeof(float),(Nr+1)*(Nth+1)*(Nph+1),outfile);
    fclose(outfile);

    /*
    strcpy(fname,"bb_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(bb_ijk,sizeof(float),(Nr+1)*(Nth+1)*(Nph+1),outfile);
    fclose(outfile);
    */

    strcpy(fname,"ta_0000.bin");
    fname[3]=48+ik;
    fname[4]=48+ih;
    fname[5]=48+id;
    fname[6]=48+iu;
    outfile = fopen(fname, "wb");
    fwrite(tau_ijk,sizeof(float),(Nr+1)*(Nth+2)*(Nph+1),outfile);
    fclose(outfile);
  }
  return(0);
}
