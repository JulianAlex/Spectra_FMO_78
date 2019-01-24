#include "spek_head.h"


void dynTheoryArray(int nz, int anz, float ***sideband, float ***zeroline,
		    fftw_complex *gt, double **gam, double *tau, double *w_m0_sh,
		    double *spek_ab, double *alpha, 
		    double *spek_sb, double **arrayspek_sb,  
		    double *spek_zl, double **arrayspek_zl,  
		    double *spek_cd, double *rotst, double *spek_ld, double *ldcff, 
		    double *spek_cdmk, double *spek_ldmk ){


  double wstep = WSTEP/WZ2WFS;  //Schrittweite w-Achse [cm-1]  
  double lineshape=0.0, linesh_sb=0.0, linesh_zl=0.0, hexp=0.0, wmo=0.0, ww=0.0;
  int nrow=0, j=0, k=0, m=0, n=0; 
  int ib[nz], ie[nz], msel[nz], nsel[nz]; 
  
  for(j=0; j<nz; j++)  //init.
  { 
    spek_ab[j]=0; spek_cd[j]=0; spek_ld[j]=0; spek_ldmk[j]=0; spek_ldmk[j]=0; 
    ib[j]=0; ie[j]=0;
    msel[j]=0; nsel[j]=0;
 }

  double **sinexct_ab = calloc( nz, sizeof(double*) );
  double **sinexct_sb = calloc( nz, sizeof(double*) );
  double **sinexct_zl = calloc( nz, sizeof(double*) );
  double **sinexct_cd = calloc( nz, sizeof(double*) );
  double **sinexct_ld = calloc( nz, sizeof(double*) );
  double **sinexct_cdmk = calloc( nz, sizeof(double*) );
  double **sinexct_ldmk = calloc( nz, sizeof(double*) );
  nrow = nz;
  while(nrow--)
  {
    sinexct_ab[nrow] = calloc( (anz+NT/2), sizeof(double) ); 
    sinexct_sb[nrow] = calloc( (anz+NT/2), sizeof(double) );
    sinexct_zl[nrow] = calloc( (anz+NT/2), sizeof(double) );
    sinexct_cd[nrow] = calloc( (anz+NT/2), sizeof(double) );
    sinexct_ld[nrow] = calloc( (anz+NT/2), sizeof(double) );
    sinexct_cdmk[nrow] = calloc( (anz+NT/2), sizeof(double) );
    sinexct_ldmk[nrow] = calloc( (anz+NT/2), sizeof(double) );
  }

  double **arrayspek_ab = calloc( nz, sizeof(double*) );
  double **arrayspek_cd = calloc( nz, sizeof(double*) );
  double **arrayspek_ld = calloc( nz, sizeof(double*) );
  double **arrayspek_cdmk = calloc( nz, sizeof(double*) );
  double **arrayspek_ldmk = calloc( nz, sizeof(double*) );
  nrow = nz;
  while(nrow--)
  {
    arrayspek_ab[nrow]   = calloc( 2*anz, sizeof(double) ); 
    arrayspek_cd[nrow]   = calloc( 2*anz, sizeof(double) ); 
    arrayspek_ld[nrow]   = calloc( 2*anz, sizeof(double) ); 
    arrayspek_cdmk[nrow] = calloc( 2*anz, sizeof(double) ); 
    arrayspek_ldmk[nrow] = calloc( 2*anz, sizeof(double) ); 
  }

  // da arrayspek_sb und _zl in main allociert werden, 
  // muessen sie hier zurueckgesetzt werden:
  for(j=0; j<nz; j++)  
    for(m=0; m < (2*anz); m++) 
    {
      arrayspek_sb[j][m] = 0.;
      arrayspek_zl[j][m] = 0.;
    }


  for(m=0; m<nz; m++)  // Index Verschiebereien
  {
    ib[m] = (int) ((SPEKMIN - w_m0_sh[m])/wstep - 0.5); //new
    ie[m] = (int) (SPEKMAX - w_m0_sh[m])/wstep;

    if(ib[m] > 0)
    {
      printf("Achtung: ib > 0!!  Something wrong with SPEKMIN SPEKMAX (head.h) %3d  %3d\n", m, ib[m]);
      printf("Exit in spek_dynTheory.c\n");
      exit(8);
    } 
    //printf("m, w_m0, ib, ie, ie-ib: %3d %7.1lf %3d %3d %3d\n", m,  w_m0_sh[m], ib[m], ie[m], ie[m]-ib[m]);
  }
   
 
  for(j=0; j<nz; j++)  // Zuordnung der passenden lineshapes, [tau]=1/cm, [DTAU]=1/fs
  {
    msel[j] = (int) (0.5 + (gam[j][j]-DGAM_0)/DGAM);    // msel[j] = (int) ( gam[j][j]/DGAM );  
    nsel[j] = (int) (0.5 + (tau[j]-DTAU_0)/DTAU);       // nsel[j] = (int) ( tau[j]/DTAU );           
      
    //printf(" %3d  %3d  %lf  %lf\n", msel[j], nsel[j], gam[j][j], tau[j]);

    if( (msel[j] < 0) || (msel[j] >= NGAM) )
    {
      printf("Achtung: Arraygrenzen!!  %3d  %3d  %lf\n", j, msel[j], gam[j][j]);
      printf("Exit in spek_dynTheory.c\n");
      exit(8);
    }

    if( nsel[j] < 0 )
    {
      printf("Achtung: DTAU_0 zu gross!!  %3d  %3d\n", j, nsel[j]);
      printf("Exit in spek_dynTheory.c\n");
      exit(8);	
    }
    if( nsel[j] >= NTAU )
    {
      if( ABS(tau[j]-TAU_0) < 1.0 )
	nsel[j] = NTAU-1; 
      else
      {
	printf("Achtung: NTAU zu klein!!  %3d  %3d\n", j, nsel[j]);
	printf("Exit in spek_dynTheory.c\n");
	exit(8);	
      }
    }
  }


  for(m=0; m<nz; m++)
  {
    wmo = w_m0_sh[m];    // [cm-1]

    hexp = exp(-creal(gt[0])*gam[m][m]);   // =Exp[-G_M(0)]

    for(k=ib[m]; k<ie[m]; k++)
    {
      ww = wmo + wstep*k;  // [cm-1]

      lineshape = hexp*( zeroline[msel[m]][nsel[m]][k+NT/2] + sideband[msel[m]][nsel[m]][k+NT/2] );  // [fs] 

      // Fuer Vergleich mit Markov, ZL(w) NICHT mit exp(-G_M(0)) multiplizieren: 

      //linesh_sb = sideband[msel[m]][nsel[m]][k+NT/2];   
      //linesh_zl = zeroline[msel[m]][nsel[m]][k+NT/2];  

      linesh_sb = hexp*sideband[msel[m]][nsel[m]][k+NT/2];   
      linesh_zl = hexp*zeroline[msel[m]][nsel[m]][k+NT/2];  

      sinexct_ab[m][k+NT/2] = ww*alpha[m]/(3*HBARC*EVWZ)*lineshape*ALPHFAC; // [e^2 Ang s]
      sinexct_sb[m][k+NT/2] = ww*alpha[m]/(3*HBARC*EVWZ)*linesh_sb*ALPHFAC; 
      sinexct_zl[m][k+NT/2] = ww*alpha[m]/(3*HBARC*EVWZ)*linesh_zl*ALPHFAC; 

      sinexct_cd[m][k+NT/2] = ww/CVAK*rotst[m]*lineshape; 
      sinexct_ld[m][k+NT/2] = ww/CVAK*ldcff[m]*lineshape; 

      sinexct_cdmk[m][k+NT/2] = ww/CVAK*rotst[m]*linesh_zl; 
      sinexct_ldmk[m][k+NT/2] = ww/CVAK*ldcff[m]*linesh_zl; 

    }
    // Konversion [Debye^2*s/m]*4.34e-12 = e^2*Ang*s, 4.34e-12=10^-10/4.8^2

    //printf("M, G_M(0)  %3d  %lf  %lf\n", m, creal(gt[0]), hexp );

  } //end m-loop

  // EINHEITEN: 
  // [alpha] = DEBYE^2, DEBYE = 4.8*e*Ang = 4.8*10^-10*e = DEBFAC*ECHARGE*ANG

  for(m=0; m<nz; m++)
    {
      n = 0;  // init.
      wmo = w_m0_sh[m];
      
      for(k = ib[m]; (wmo + wstep*(k+1)) < SPEKMIN; k++ );
      
      for(k = k; (wmo + wstep*(k-1)) < SPEKMAX; k++ )
	{
	  arrayspek_ab[m][n]   = sinexct_ab[m][k+NT/2]; 
	  arrayspek_sb[m][n]   = sinexct_sb[m][k+NT/2];
	  arrayspek_zl[m][n]   = sinexct_zl[m][k+NT/2]; 
	  arrayspek_cd[m][n]   = sinexct_cd[m][k+NT/2];
	  arrayspek_ld[m][n]   = sinexct_ld[m][k+NT/2];
	  arrayspek_cdmk[m][n] = sinexct_cdmk[m][k+NT/2];
	  arrayspek_ldmk[m][n] = sinexct_ldmk[m][k+NT/2];
 	  n++;
	}
    } 

  
  //  for(m=0; m<nz; m++)
  //    for(j = 0; j<anz; j++)
  //      printf("%3d %3d %12.6lf %12.6lf %12.6lf \n", m, j, arrayspek_ab[m][j], arrayspek_sb[m][j], arrayspek_zl[m][j]);
  


  // Spektren fuer alle Pigmente zusammen: 

  for(j = 0; j<anz; j++)
  {
    spek_ab[j] = 0.0; spek_sb[j] = 0.0; spek_zl[j] = 0.0; spek_cd[j] = 0.0; spek_ld[j] = 0.0; 
    spek_cdmk[j] = 0.0; spek_ldmk[j] = 0.0;

    for(m=0; m<nz; m++)
    {
      spek_ab[j]   += arrayspek_ab[m][j];
      spek_sb[j]   += arrayspek_sb[m][j];
      spek_zl[j]   += arrayspek_zl[m][j];
      spek_cd[j]   += arrayspek_cd[m][j];
      spek_ld[j]   += arrayspek_ld[m][j];
      spek_cdmk[j] += arrayspek_cdmk[m][j];
      spek_ldmk[j] += arrayspek_ldmk[m][j];
    }
  }

  //-------------------------------------


  nrow = nz;
  while(nrow--)
  {
    free(arrayspek_ab[nrow]);
    free(arrayspek_cd[nrow]);
    free(arrayspek_ld[nrow]);
    free(arrayspek_cdmk[nrow]);
    free(arrayspek_ldmk[nrow]);
  }
  free(arrayspek_ab); arrayspek_ab = NULL;
  free(arrayspek_cd); arrayspek_cd = NULL;
  free(arrayspek_ld); arrayspek_ld = NULL;
  free(arrayspek_cdmk); arrayspek_cdmk = NULL;
  free(arrayspek_ldmk); arrayspek_ldmk = NULL;

  nrow = nz;
  while(nrow--)
  {
    free(sinexct_ab[nrow]);
    free(sinexct_sb[nrow]);
    free(sinexct_zl[nrow]);
    free(sinexct_cd[nrow]);
    free(sinexct_ld[nrow]);
    free(sinexct_cdmk[nrow]);
    free(sinexct_ldmk[nrow]);
  }
  free(sinexct_ab); sinexct_ab = NULL;
  free(sinexct_sb); sinexct_sb = NULL;
  free(sinexct_zl); sinexct_zl = NULL;
  free(sinexct_cd); sinexct_cd = NULL;
  free(sinexct_ld); sinexct_ld = NULL;
  free(sinexct_cdmk); sinexct_cdmk = NULL;
  free(sinexct_ldmk); sinexct_ldmk = NULL;

}//end


//================================================================================


void markovApprox(int nz, int anz, double *alpha, double *rotst, double *ld_koeff, 
		  double *tau, double *w_m0_sh, 
		  double *nmark_od, double *nmark_cd, double *nmark_ld){
 
  // Bem: Die Markov-Näherung ist auch eine dynamische Theorie

  int j, k;

  double w, step = WSTEP/WZ2WFS;
  double f[nz][anz], g[nz][anz], h[nz][anz];

  for(k=0; k<nz; k++){  
    for(j=0; j<anz; j++){ 
      w = SPEKMIN + j*step;   // w in wellenzahlen !!

      f[k][j] =    w*pow(10,-8)*alpha[k]/tau[k]/( SQ((w - w_m0_sh[k])*WZ2WFS) + SQ(1/tau[k]) );
      g[k][j] =    w*pow(10,-8)*rotst[k]/tau[k]/( SQ((w - w_m0_sh[k])*WZ2WFS) + SQ(1/tau[k]) );
      h[k][j] = w*pow(10,-8)*ld_koeff[k]/tau[k]/( SQ((w - w_m0_sh[k])*WZ2WFS) + SQ(1/tau[k]) );

    }
  } 



  for(j=0; j<anz; j++)
  {  
    for(k=0; k<nz; k++)
    {  
      nmark_od[j] += f[k][j]; 
      nmark_cd[j] += g[k][j];
      nmark_ld[j] += h[k][j];
    }          
  }

 
}


//===================================================================================


void dynTheoryArrayHB(int nz, int anz, float ***sideband, float ***zeroline,
		      fftw_complex *gt, double **gam, double *tau, double *w_m0_sh,
		      double *spek_hb, double *alpha, double **vab, double **vab_hb){

  //  FILE *out_file; char *out_dat;


  double wstep = WSTEP/WZ2WFS;  //Schrittweite w-Achse [cm-1]  
  double lineshape=0.0, hexp=0.0, wmo=0.0, ww=0.0;
  int nrow=0, j=0, k=0, m=0, n=0; 
  int ib[nz], ie[nz], msel[nz], nsel[nz]; 
  
  for(j=0; j<nz; j++)  //init.
  { 
    spek_hb[j]=0; 
    ib[j]=0; ie[j]=0;
    msel[j]=0; nsel[j]=0;
  }

  double **sinexct_hb   = calloc( nz, sizeof(double*) );
  double **arrayspek_hb = calloc( nz, sizeof(double*) );
  nrow = nz;
  while(nrow--)
  {
    sinexct_hb[nrow]   = calloc( (anz+NT/2+1), sizeof(double) ); 
    arrayspek_hb[nrow] = calloc( 4*anz, sizeof(double) ); 
  }


  for(m=0; m<nz; m++)  // Index Verschiebereien
  {
    // ib[m] = (int) (SPEKMIN - w_m0_sh[m])/wstep; // wrong

    ib[m] = (int) ((SPEKMIN - w_m0_sh[m])/wstep - 0.5); //new
    ie[m] = (int) (SPEKMAX - w_m0_sh[m])/wstep;

    //printf("m, w_m0, ib, ie, ie-ib: %3d %7.0lf %3d %3d %3d\n", m,  w_m0_sh[m], ib[m], ie[m], ie[m]-ib[m]);
  }
   
 
  for(j=0; j<nz; j++)  // Zuordnung der passenden lineshapes, [tau]=1/cm, [DTAU]=1/fs
  {
    msel[j] = (int) (0.5+(gam[j][j]-DGAM_0)/DGAM);    // msel[j] = (int) ( gam[j][j]/DGAM );  
    nsel[j] = (int) (0.5+(tau[j]-DTAU_0)/DTAU);       // nsel[j] = (int) ( tau[j]/DTAU );           
      
    //printf(" %3d  %3d  %lf  %lf\n", msel[j], nsel[j], gam[j][j], tau[j]);

    if( (msel[j] < 0) || (msel[j] >= NGAM) )
    {
      printf("Achtung: Arraygrenzen!!  %3d  %3d  %lf\n", j, msel[j], gam[j][j]);
      printf("Exit in spek_dynTheory.c\n");
      exit(8);
    }

    if( nsel[j] < 0 )
    {
      printf("Achtung: DTAU_0 zu gross!!  %3d  %3d\n", j, nsel[j]);
      printf("Exit in spek_dynTheory.c\n");
      exit(8);	
    }
    if( nsel[j] >= NTAU )
    {
      if( ABS(tau[j]-TAU_0) < 1.0 )
	nsel[j] = NTAU-1; 
      else
      {
	printf("Achtung: NTAU zu klein!!  %3d  %3d\n", j, nsel[j]);
	printf("Exit in spek_dynTheory.c\n");
	exit(8);	
      }
    }
  }


  for(m=0; m<nz; m++)
  {
    wmo = w_m0_sh[m];    // [cm-1]

    hexp = exp(-creal(gt[0])*gam[m][m]);   // =Exp[-G_M(0)]

    //del_e = ABS(E_BURN - vab_hb[m][m]);  //ABS(vab[m][m]-vab_hb[m][m]);
    //belohnfkt = (1.0 + a*exp(-b*del_e)); 
    //printf("%6.0lf %6.0lf %5.2lf %g\n", vab[m][m], vab_hb[m][m], del_e, belohnfkt);

    for(k=ib[m]; k<ie[m]; k++)
    {
      ww = wmo + wstep*k;  // [cm-1]

      lineshape = hexp*( zeroline[msel[m]][nsel[m]][k+NT/2] + sideband[msel[m]][nsel[m]][k+NT/2] );  // [fs] 

      sinexct_hb[m][k+NT/2] = 4*PISQ*ww*NINDEX*alpha[m]/(3*HBARC*EVWZ)*lineshape*ALPHFAC; // [e^2 Ang s]
    }
    // Konversion [Debye^2*s/m]*4.34e-12 = e^2*Ang*s, 4.34e-12=10^-10/4.8^2

    //printf("M, G_M(0)  %3d  %lf  %lf\n", m, creal(gt[0]), hexp );

  } //end m-loop

  // EINHEITEN: 
  // [alpha] = DEBYE^2, DEBYE = 4.8*e*Ang = 4.8*10^-10*e = DEBFAC*ECHARGE*ANG

  for(m=0; m<nz; m++)
    {
      n = 0;  // init.
      wmo = w_m0_sh[m];
      
      for(k = ib[m]; (wmo + wstep*(k+1)) < SPEKMIN; k++ );
      
      for(k = k; (wmo + wstep*(k-1)) < SPEKMAX; k++ )
	{
	  arrayspek_hb[m][n] = sinexct_hb[m][k+NT/2]; 
	  //printf("%3d %3d %3d %12.6lf %12.6lf\n", m, n, k, arrayspek_hb[m][n], sinexct_hb[m][k]);
	  n++;
	}
    } 


  // Spektren fuer alle Exzitonen zusammen: 

  for(j = 0; j<anz; j++)
  {
    spek_hb[j] = 0.0;
    for(m=0; m<nz; m++)
    {
      spek_hb[j] += arrayspek_hb[m][j];
    }
  }

  //-------------------------------------

  nrow = nz;
  while(nrow--)
  {
    free(sinexct_hb[nrow]);
    free(arrayspek_hb[nrow]);
  }
  free(sinexct_hb); sinexct_hb = NULL;
  free(arrayspek_hb); arrayspek_hb = NULL;


}//end dyntheo.c

















/********************************************************************

wie testen ob  exp(G_M(0)) - 1 = int(SB)/int(ZL)  ?

T = 0.01
infile_*.dat alle bis auf ein pigment loeschen (nz = 1)
out_nonMarkov.dat nur fuer zeroline und nur fuer sideband vergleichen
normierung ausschalten!

*********************************************************************/




  /*
  out_dat="out_lineshifted.dat";
  out_file = fopen(out_dat, "w");

  for(m=1; m<=nz; m++)
    for(k=ib[m]; k<ie[m]; k++){
      ww = w_m0_sh[m-1]*WZ2WFS+TWOPI*k/(NT*DELTA); //in 1/fs
      fprintf(out_file," %d  %lf  %lf\n", k, ww/WZ2WFS, sinexct_ab[m-1][k]);
    }
  fclose(out_file);




  for(j=0; j<nz; j++) 
    for(m=0; m < (2*anz); m++) 
    {
      arrayspek_sb[j][m] = 0.;
      arrayspek_zl[j][m] = 0.;
    }



  */
