#include "spek_head.h"


//===========================================================================

int outputStickSpek(double *nr_ew, double *alpha, int nz){

  FILE *out_file;
  char *out_dat="out_stick.dat";
  int j=0;


  out_file = fopen(out_dat, "w");
  if((out_file = fopen(out_dat, "w")) == NULL){
    printf( "fopen error\n" );
    return 1;
  }

  for(j=0; j<nz; j++)
    fprintf(out_file,"%lf    %lf \n",
	    nr_ew[j+1], alpha[j]/10);
  fclose(out_file);

  return 0;

}


//===========================================================================

int outputSpectrum(double *sumspek, int anz, char *out_dat){


  FILE *out_file;
  int k=0; 
  double w=0.0, wstep=0.0;

  wstep = TWOPI/(NT*DELTA*WZ2WFS); //[cm-1] 

  out_file = fopen(out_dat, "w");

  if((out_file = fopen(out_dat, "w")) == NULL)
  {
    printf( "%s fopen error\n", out_dat);
    return 1;
  }

  for(k=0; k<anz; k++)
  {             
    w = SPEKMIN + k*wstep;                                        //- wstep/2 + k*wstep;
    fprintf(out_file, "%lf  %lf  %lf \n", w, 1e+7/w, sumspek[k]);  // [cm-1], [nm]
  }
  fclose(out_file);

  return 0;

}
//===========================================================================

int outputSpectrumShift(double *sumspek, char *out_dat){


  FILE *out_file;
  int k=0, anz=0;
  double w=0.0, wstep=0.0;

  wstep = TWOPI/(NT*DELTA*WZ2WFS); //[cm-1]
  anz   = (int) (SPEKDIFF/wstep + 1);    
  printf("Anz %d  %lf\n", anz, wstep);

  out_file = fopen(out_dat, "w");

  if((out_file = fopen(out_dat, "w")) == NULL)
  {
    printf( "%s fopen error\n", out_dat);
    return 1;
  }

  for(k=0; k<=anz+1; k++)
  {             
    w = SPEKMIN  + k*wstep + 20;                                        //- wstep/2 + k*wstep;
    fprintf(out_file, "%lf  %lf  %lf \n", w, 1e+7/w, sumspek[k]);  // [cm-1], [nm]
  }
  fclose(out_file);

  return 0;

}

//=======================================================================

int outputPigmExc(double **npig, double **nexc){

  FILE *out_file;
  char *out_dat="out_pigment.dat";
  int j=0;
  double w=0.0;

  out_file = fopen(out_dat, "w");
  if((out_file = fopen(out_dat, "w")) == NULL){
    printf( "%s fopen error\n", out_dat);
    return 1;
  }
  for(j=0; j<NEXC; j++){
    w = SPEKMIN+j*EXCLENG;
    //    fprintf(out_file,"%lf   %g   %g   %g   %g   %g   %g   %g\n", w,
    //	    npig[0][j], npig[1][j], npig[2][j], npig[3][j], npig[4][j], npig[5][j], npig[6][j]);
    fprintf(out_file," %lf   %g  %g \n", w,
    	    npig[0][j], npig[1][j]);
  }
  fclose(out_file);

  out_dat="out_exciton.dat";

  out_file = fopen(out_dat, "w");
  if((out_file = fopen(out_dat, "w")) == NULL){
    printf( "%s fopen error\n", out_dat);
    return 1;
  }
  for(j=0; j<NEXC; j++){
    w = SPEKMIN+j*EXCLENG;
    fprintf(out_file,"%lf   %lf   %lf \n", w,
	    nexc[0][j], nexc[1][j]);
  }
  fclose(out_file);

  return 0;

}



