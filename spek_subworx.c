// Julian Adolphs
// File "koorParaInter.c" enhaelt die Proceduren:

// void normierung
// void normSB
// void sumMC
// void koord
// void dipolmom
// void pigmCenters
// void parameter
// void interaction
// void calcDipoleStrength
// void calcDipoleMoment
// void coulombCouplings
// void matrix2dia
// void ewev

#include "spek_head.h"

// ==========================================================================

void normierung(int n, double leng, double *ab){

  int j=0;
  double sum_ab=0.0;


  // FlaechenNormierung

  for(j=0; j<n; j++){
    sum_ab += ABS(ab[j]);
  }

  for(j=0; j<n; j++){
    ab[j] *= 100/(sum_ab*leng);
    //printf("%d %8.6lf\n", j, ab[j]);
  }

}

//============================================================================

void normSB(double *ab,double *hb)
{

  int j=0;
  double sum_ab=0.0;

  double wstep = TWOPI/(NT*DELTA*WZ2WFS);
  int    anz   = (int) (SPEKDIFF/wstep + 1);    

  // FlaechenNormierung

  for(j=0; j < anz; j++)
  {
    sum_ab   += ABS(ab[j]);
  }

  for(j=0; j < anz; j++)
  {
    if (sum_ab > 0.0)
        ab[j]  *= 100.0/(sum_ab*wstep);
        hb[j]  *= 100.0/(sum_ab*wstep);
    //printf("%d %8.6lf\n", j, ab[j]);
  }

}

//=====================================================================================

void koord(int nz, double *x, double *y, double *z, 
	   double **na, double **nb, double **nc, double **nd){
  int i=0, j=0;

  for(i=0;i<nz;i++)
    { 
      j=i*NNATM;
      na[i][0]=x[j];   na[i][1]=y[j];   na[i][2]=z[j];
      nb[i][0]=x[j+1]; nb[i][1]=y[j+1]; nb[i][2]=z[j+1];
      nc[i][0]=x[j+2]; nc[i][1]=y[j+2]; nc[i][2]=z[j+2];
      nd[i][0]=x[j+3]; nd[i][1]=y[j+3]; nd[i][2]=z[j+3];
    }
}

//============================================================================

 void trimSymAx(int nz, double **trim_r, double *sym_ax){

   // Funktion "trim_sym_ax.c
   // Berechnet Symmetrie-Axe des Trimers, notwendig fuer LD-Spektrum

  int i;
  double abs=0.0, p1[3], p2[3], p3[3], a[3], b[3], n[3];

  // waehle zum auffinden der symmetrieachse die koordinaten von Mg in BCL 1
  // jeweils in monomer A, B und C => Punkte p1, p2, p3. 
  // SymmAxe ist dann senkrecht auf dieser Ebene.
  // vec_n = vec_p_12 x vec_p_13 ,  vec_p_12 = p_2 - p_1  

  for(i=0; i<3; i++){
    p1[i] = trim_r[0][i];
    p2[i] = trim_r[nz][i];
    p3[i] = trim_r[2*nz][i];
 
    a[i] = p2[i] - p1[i];
    b[i] = p3[i] - p1[i];

    //printf(" %f  %f  %f\n", p1[i], p2[i], p3[i]);
  }

  // Vektorprodukt:

  n[0] = a[1]*b[2] - a[2]*b[1];
  n[1] = a[2]*b[0] - a[0]*b[2];
  n[2] = a[0]*b[1] - a[1]*b[0];

  // Normierung:

  for(i=0; i<3 ; i++){
    abs += SQ(n[i]);
  }  
  abs = sqrt(abs);

  for(i=0; i<3 ; i++)
    sym_ax[i] = n[i]/abs; 

}

//=====================================================================================

// Berechnet normierte Richtungvektoren der Dipolmoment qx = NA-NC und qy=NB-ND 
// Einheit: [qy] = 1

void dipolmom(int nz, double **na, double **nb, double **nc,
	      double **nd, double **qx, double **qy){

  int i=0, j=0;
  double abs_x=0.0, abs_y=0.0, beta = BETA*PI/180.0;

  for(i=0;i<nz;i++){
    for(j=0;j<3;j++){
      qx[i][j] = na[i][j]-nc[i][j];
      qy[i][j] = nb[i][j]-nd[i][j];
    }
  }

  // Normieren mit Betrag!

  for(i=0;i<nz;i++){
    abs_x=0; abs_y=0;
    for(j=0;j<3;j++){
      abs_x += SQ(qx[i][j]);
      abs_y += SQ(qy[i][j]);
    }
    for(j=0;j<3;j++){
      qx[i][j] = qx[i][j]/sqrt(abs_x);
      qy[i][j] = qy[i][j]/sqrt(abs_y);
    }
  }

  // Drehung um Winkel Beta

  for(i=0; i<nz; i++)
    for(j=0; j<3; j++)
      qy[i][j] = qx[i][j]*sin(beta) + qy[i][j]*cos(beta);

  printf("\n qy \n");

  for(i=0;i<nz;i++){
    for(j=0;j<3;j++)
      printf(" %9.6f ", qy[i][j]);
    printf(" \n");
  }
  printf(" \n");
  
}

//============================================================================

void pigmCenters(int nz, double **mg, double **nb, double **nd,
		double **delta, double **rij){

  // Bestimmung der chl-Zentren,
  // Abstandsbetrag rij[i][j] zwischen i und j,
  // normierter Richtungsvektor delta[i][j][k], i,j=nz, k=3

  int i=0, j=0, k=0;
  double abs=0.0, vektor[3];

  for(j=0; j<3; j++) vektor[j]=0.0;


  for(i=0; i<nz; i++){
    for(j=0; j<3; j++){
      vektor[j]= nd[i][j] - nb[i][j];
      mg[i][j] = nb[i][j] + vektor[j]/2.0; 
    }
  }


  // Berechne Differenzvektoren zwischen chl-Zentren
  // Delta eigentlich 3dim Matrix delta[i][j][k] mit i,j=0,...,nz-1, k=0,1,2
  // Mit der Ersetzungsvorschrift delta[i][j][k] => delta[i+j*nz][k]

  for(i=0; i<nz-1; i++)
    {
      for(j=i+1; j<nz; j++)
	{
	  for(k=0; k<3; k++)
	    delta[i+j*nz][k] = mg[i][k]-mg[j][k];
	}
    }
  
  //  rij ist Abstands-Betrag (in Angstr) zwischen Mg_i und Mg_j

    for(i=0; i<nz-1; i++){
      for(j=i+1; j<nz; j++){
	abs = 0;
	for(k=0; k<3; k++)
	  abs += SQ(delta[i+j*nz][k]);
	abs = sqrt(abs);
	rij[i][j] = abs;
	rij[j][i] = rij[i][j];
	for(k=0; k<3; k++)
	  delta[i+j*nz][k] = delta[i+j*nz][k]/abs;
      }
      rij[i][i]=0.0;
    }

    for(i=0; i<nz-1; i++){
      for(j=i+1; j<nz; j++){
	for(k=0; k<3; k++)
	  delta[j+i*nz][k] = -delta[i+j*nz][k] ;
      }
    }

}

//============================================================================

void parameter(double *site, double *sigma, double *dipst, int bclnum, int nz){

  int i=0;

  // site energies in wave-numbers (1/cm)

  site[0] = SITE_0;
  site[1] = SITE_1;
  site[2] = SITE_2;
  site[3] = SITE_3;
  site[4] = SITE_4;
  site[5] = SITE_5;
  site[6] = SITE_6;
  site[7] = SITE_7;

  // Width of Gauss-Distrib of site-Energies in MonteCarlo

  sigma[0] = 0.425*FWHM_0;   // SIGMA = 0.425 FWHM
  sigma[1] = 0.425*FWHM_1;
  sigma[2] = 0.425*FWHM_2;
  sigma[3] = 0.425*FWHM_3;
  sigma[4] = 0.425*FWHM_4;
  sigma[5] = 0.425*FWHM_5;
  sigma[6] = 0.425*FWHM_6;
  sigma[7] = 0.425*FWHM_7;
 
  for(i=0; i<nz; i++){ 
    dipst[i] = DIPST;
  }

  if(bclnum == 8)
  {
     dipst[7] = DIPST;
     printf("FMO variant with 8 BCHls\n\n"); 
  }
  else if(bclnum == 7)
  {     
    dipst[7] = 0.0;
    printf("FMO variant with 7 BCHls\n\n");  
  }
  else
  {
    printf("BCl-Number wrong! %d n\n", bclnum); 
    exit(8);
  }

  printf(" \n");
  for(i=0; i<nz; i++)
    printf(" %6.0lf", site[i]); printf(" \n");
  printf(" \n");

  for(i=0; i<nz; i++)
    printf(" %6.0lf", sigma[i]/0.425);
  printf(" \n"); printf(" \n");
}

//============================================================================

void interaction(int nz, double **qy, double **delta, double **vab,
		 double **rij, double *dipst){

  // Berechnet Wechselwirkungsenergie zweier Dipole

  int i=0, j=1, k=0;
  double h0=0.0, h1=0.0, h2=0.0, hvek[nz][nz][3];

  for(i=0; i<nz; i++) 
    for(j=0; j<nz; j++) 
      for(k=0; k<3; k++) 
	hvek[i][j][k]=0.0;

  printf("Dipole strength: %6.3lf,   effektive dipole strength %6.3lf\n\n",
	 dipst[0], sqrt(DIELEK)*dipst[0]);

  for(i=0; i<nz-1; i++){
    for(j=i+1; j<nz; j++){
      h0=0; h1=0; h2=0;
      for(k=0; k<3; k++){
	h0 += delta[i+j*nz][k]*qy[i][k];
        h1 += delta[i+j*nz][k]*qy[j][k];
	h2 += qy[i][k]*qy[j][k];
      }
      hvek[i][j][0]=h0;
      hvek[i][j][1]=h1;
      hvek[i][j][2]=h2;
    }
  }

  // Umrechnung der Energie(D^2/A^3) in E(1/cm) mit Faktor 5040

  for(i=0; i<nz-1; i++)
  {
    for(j=i+1; j<nz; j++)
    {
      vab[i][j] = hvek[i][j][2]-3*hvek[i][j][1]*hvek[i][j][0];
      vab[i][j] *= DIELEK/(pow((rij[i][j]), 3))*5040.84*dipst[i]*dipst[j];
    }
  }

  for(i=0; i<nz; i++)
    vab[i][i] = 0.0;

  for(i=1; i<nz; i++){
    for(j=0; j<i; j++){
      vab[i][j] = vab[j][i];
    }
  }

  for(i=0; i<nz; i++)
  {
    for(j=0; j<i; j++)
    {
      printf("%d %d %11.4f \n", i, j, vab[i][j]);
    }
    printf("\n");
  }
  printf("\n");

}

//============================================================================

void calcDipoleStrength(int n1, double *x1, double *y1, double *z1,
			double *charge){

  int i=0;
  double dx[n1], dy[n1], dz[n1], corrfac=0.0;

  for(i=0; i<n1; i++){ dx[i]=0.0; dy[i]=0.0; dz[i]=0.0; }


  // Dipolstaerke: e*Ang = 4.8 Debye

  //  for(i=0; i<n1; i++)
  //    charge[i] *= MILLI;

  double qx=0.0, qy=0.0, qz=0.0, dipst=0.0;

  for(i=0; i<n1; i++){
    dx[i] = DEBFAC*charge[i]*x1[i];
    dy[i] = DEBFAC*charge[i]*y1[i];
    dz[i] = DEBFAC*charge[i]*z1[i];
  }

  for(i=0; i<n1; i++){
    qx += dx[i];
    qy += dy[i];
    qz += dz[i];
  }

  dipst = sqrt(SQ(qx)+SQ(qy)+SQ(qz));

  printf("%10.6lf  dipst/Debey\n", dipst);
  printf("\n");

  corrfac = DIPST/dipst;

  // === Probe ====

  for(i=0; i<n1; i++){
    dx[i] *= corrfac;
    dy[i] *= corrfac;
    dz[i] *= corrfac;
  }

  qx = 0.0; qy = 0.0; qz = 0.0;
  for(i=0; i<n1; i++){
    qx += dx[i];
    qy += dy[i];
    qz += dz[i];
  }

  dipst = sqrt(SQ(qx)+SQ(qy)+SQ(qz));

  printf("%10.6lf  renorm.dipst/Debey\n", dipst);
  printf("\n");
  printf("%10.6lf  Renormierungs-Faktor\n", corrfac);
  printf("\n");

  //  for(i=0; i<n1; i++)
  //    printf("%d %lf\n", i, charge[i]);
  //  printf("\n");

  for(i=0; i<n1; i++)
    charge[i] *= corrfac;

  //  for(i=0; i<n1; i++)
  //    printf("%d %lf\n", i, charge[i]);
  //  printf("\n");

}

//=================================================================

void calcDipoleMoment(int nz, int n1, double *x7, double *y7, double *z7,
		      double *charge, double **qy){

  int i=0, j=0, k=0;
  double dx[n1], dy[n1], dz[n1];

  for(i=0; i<nz; i++){ dx[i]=0.0; dy[i]=0.0; dz[i]=0.0; }

  // Dipolstaerke: e*Ang = 4.8 Debye

  double qx1=0.0, qx2=0.0, qx3=0.0, dipst=0.0;

  for(i=0; i<nz; i++){
    for( j = k*n1; j < (k+1)*n1; j++ ){
      dx[j%n1] = DEBFAC*charge[j%n1]*x7[j];
      dy[j%n1] = DEBFAC*charge[j%n1]*y7[j];
      dz[j%n1] = DEBFAC*charge[j%n1]*z7[j];
    }
    k++;
    qx1 = 0.0; qx2 = 0.0; qx3 = 0.0;
    for(j=0; j<n1; j++){
      qx1 += dx[j];
      qx2 += dy[j];
      qx3 += dz[j];
    }
    dipst = sqrt(SQ(qx1)+SQ(qx2)+SQ(qx3));
    qy[i][0] = qx1/dipst;
    qy[i][1] = qx2/dipst;
    qy[i][2] = qx3/dipst;
  }

  printf("Dipole-Moment:\n\n");
  for(i=0; i<nz; i++)
    printf("%d  %10.6lf  %10.6lf  %10.6lf\n", i, qy[i][0], qy[i][1], qy[i][2]);
  printf("\n");

}

//=======================================================================

void coulombCouplings(int nz, int n1, double *x7, double *y7, double *z7,
		      double *charge, double **vab){

  // Ausrechnen der Coulomb-Kopplungen aus transCharges

  int i=0, j=0, k=0, m=0;
  double dx=0.0, dy=0.0, dz=0.0, dr=0.0, sum=0.0;

  printf("WW-Energie im Vakuum/cm-1:\n");
  printf("\n");

  for(i=0; i<nz; i++){
    for(j=0; j<nz; j++){
      if( j != i ){
	sum = 0.0;
	for(k=0; k<n1; k++){
	  for(m=0; m<n1; m++){
	    dx = x7[i*n1+k] - x7[j*n1+m];
	    dy = y7[i*n1+k] - y7[j*n1+m];
	    dz = z7[i*n1+k] - z7[j*n1+m];
	    dr = sqrt( dx*dx + dy*dy + dz*dz );
	    sum += charge[k]*charge[m]/dr;
	  } // end-m
	} // end-k
      } // endif
      else if( j == i ){
	sum = 0.0;
      }
      vab[i][j] = sum*1000*ECHARGE*EVWZ/(4*PI*EPS_0*ANG);
      printf(" %8.3lf", vab[i][j]);
    } // end-j
    printf("\n");
  } // end-i

  printf("\n");
  printf("WW-Energie im Dielektrikum/cm-1:\n");
  printf("\n");

  for(i=0; i<nz; i++){
    for(j=0; j<nz; j++){
      vab[i][j] *= DIELEK*DIPFAK;
      printf(" %8.3lf", vab[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  printf("WW-Energie im Dielektrikum/eV:\n");
  printf("\n");

  for(i=0; i<nz; i++){
    for(j=0; j<nz; j++){
      printf(" %8.3lf", vab[i][j]/EVWZ);
    }
    printf("\n");
  }
  printf("\n");

}

//================================================================================================

void matrix2diaHB(int nz, double *site, double *sigma, int *seed, double *w_m0_sh, double **vab, double **vab_hb, 
		   double *hamilt_hb, int pigm_no)
{
  // Hier wird jeweils eine Site-Energy des gebrannten Pigments pigm_no neu gewuerfelt
   
  // FILE *out_file; char *out_dat; out_dat="out_random_site.dat"; out_file = fopen(out_dat, "a");

  int j=0, k=0;

  for(k=0; k<nz; k++)    //uebernehme urspruengliche Matrix
  {
    for(j=0; j<nz; j++)
    {
      vab_hb[k][j] = vab[k][j];
    }
  }


  // Neu Würfeln der Site Energie des Pigments pigm_no 

  vab_hb[pigm_no][pigm_no] = r8_normal_ab(site[pigm_no], sigma[pigm_no], seed); 

  for(k=0; k<nz; k++)
  {
    for(j=0; j<nz; j++)
    {
      hamilt_hb[k*nz+j] = vab_hb[k][j];
    }
  }
    
}

//================================================================================================

void matrix2diaHBres(int nz, double *site, double *sigma, int *seed, double *w_m0_sh, double **vab, double **vab_hb, 
		   double *hamilt_hb, int pigm_no)
{
  // Hier wird jeweils eine Site-Energy des gebrannten Pigments pigm_no neu gewuerfelt
   
  // FILE *out_file; char *out_dat; out_dat="out_random_site.dat"; out_file = fopen(out_dat, "a");

  int j=0, k=0;    

  for(k=0; k<nz; k++)    //uebernehme urspruengliche Matrix
  {
    for(j=0; j<nz; j++)
    {
      vab_hb[k][j] = vab[k][j];
    }
  }

  // Neu Würfeln der Site Energie des Pigments pigm_no 
  // wenn in ZL des energetisch tiefsten Exc.Zstds |M=1> gebrannt wird

  vab_hb[pigm_no][pigm_no] = r8_normal_ab( vab[pigm_no][pigm_no], WBZL, seed ); 


  for(k=0; k<nz; k++)
  {
    for(j=0; j<nz; j++)
    {
      hamilt_hb[k*nz+j] = vab_hb[k][j];
    }
  }
    
}

//====================================================================================================

void matrix2dia(int nz, double *site, double *sigma, int *seed, double **vab, double *hamilt)
{
  // WW-Matrix

  int j=0, k=0, nmc=N;
  double site_h[nz];   
  
  // FILE *out_file; char *out_dat; out_dat="out_random_site.dat"; out_file = fopen(out_dat, "a");
  
  for(j=0; j<nz; j++) 
    site_h[j]=0.0; 

  if(nmc == 1)
  {
    for(j=0; j<nz; j++){
      site_h[j] = site[j];  //wenn N=1 soll homog Spektrum berechnet werden
    }
  }
  else
  {
    for(j=0; j<nz; j++)
      {
	site_h[j] =  r8_normal_ab(site[j], sigma[j], seed); //Würfeln der Site Energies (Normalverteilung)
      }
  }

  for(j=0; j<nz; j++)
    vab[j][j] = site_h[j];  //Diag Elemente = Site Energies, Off Diag Elemente standen vorher schon drinnen (wurden eingelesen)

  for(k=0; k<nz; k++){
    for(j=0; j<nz; j++){
      hamilt[k*nz+j] = vab[k][j];
    }
  }

  //  fclose(out_file);
}

//============================================================================


void ewev(int nz, double *eigval, double *ev_h, double **eigvec){


  int j=0, k=0;

  for(j=0; j<nz; j++){
    for(k=0; k<nz; k++){
      eigvec[j][k] = ev_h[j+k*nz];  
    }
  }

  // anscheinend gibt eispack die ev zeilen- und spalten-vertauscht aus

  /***************************
  printf(" Eigenwerte: ");
  for(j=0;j<nz;j++)
    printf(" %10.2f", eigval[j]);
  printf("\n");


  printf(" Eigenvektoren: \n");
  for(j=0;j<nz;j++){
    for(k=0; k<nz; k++){
      printf(" %12.9f", eigvec[j][k]);
    }
    printf("\n");
  }
  printf("\n");
  ****************************/

}

//============================================================================

// Linearer Absorptionskoeffizient

void linAbsKoeff(int nz, double **qy, double **eigvec, double *alpha, double *dipst){

  int i=0, j=0, k=0;
  double alpha_h=0.0;

  for(k=0; k<nz; k++)
    alpha[k] = 0.0;

  for(k=0; k<nz; k++){
    for(j=0; j<3; j++){
      alpha_h = 0.0;
      for(i=0; i<nz; i++)
	alpha_h += eigvec[i][k]*qy[i][j]*dipst[i]; // berechne hier mu
      alpha[k] += SQ(alpha_h);                        // alpha = mu^2
    }
  }

  /*
  printf("alpha:     ");
  for(k=0; k<nz; k++)
    printf("  %10.2lf", alpha[k]);
  printf("\n\n");
  */

}


//============================================================================

// Rotational Strength

void rotatStrength(int nz, double **qy, double **eigvec, double *rotst, 
		   double *dipst, double **rij, double **delta){

  int i, j, k, m;  
  double vecpro[nz][nz][3], spat[nz][nz]; //vector-product

  for(k=0; k<nz; k++) 
    rotst[k] = 0.0;

  for(k=0; k<nz; k++) 
    for(m=0; m<nz; m++)
      spat[k][m] = 0.0; 

  for(i=0; i<nz; i++)
  {   
    for(j=0; j<nz; j++)
    {
      vecpro[i][j][0] = ( qy[i][1]*qy[j][2] - qy[i][2]*qy[j][1] )*dipst[i]*dipst[j]; 
      vecpro[i][j][1] = ( qy[i][2]*qy[j][0] - qy[i][0]*qy[j][2] )*dipst[i]*dipst[j];
      vecpro[i][j][2] = ( qy[i][0]*qy[j][1] - qy[i][1]*qy[j][0] )*dipst[i]*dipst[j]; 
    }
  }
    

  // (Spat-) Produkt  R_mn.(mu_m x mu_n),    (R_mn=rij*delta, mu_m x mu_n=vecpro) 

  for(i=0; i<nz; i++)
  {   
    for(j=0; j<nz; j++)
    {
      for(k=0; k<3; k++)
      {
	spat[i][j] += rij[i][j]*delta[i+j*nz][k]*vecpro[i][j][k];
      }  
      //printf("spat  %d %d %lf\n", i,j, spat[i][j]); 
    }
  }
  
  for(m=0; m<nz; m++)
  {
    for(i=0; i<nz; i++)  
    {  
      for(j=0; j<i; j++)
      {
	rotst[m] += spat[i][j]*eigvec[i][m]*eigvec[j][m];

	//printf(" %3d %3d %3d %lf\n", m, i, j, rotst[m]); 
      }
    }
  }

      
}

//============================================================================

//   LD = 0.5 * alpha * ( 1 - 3(cos(Theta_K))^2 ),   abs(mu_K)^2 = alpha
//   cos(Theta_K) = sym_ax*qy[K] 
//   LD = LD_perp - LD_para      

void linDikroKoeff(int nz, double *sym_ax, double *alpha, double **qy, 
		   double *ld_koeff, double **eigvec, double *dipst){

  int i, j, k;
  double cos_theta[nz], mu_koll[nz][3];
 
  for(k=0; k<nz; k++){  //Initialisieren
    cos_theta[k] = 0;
    for(j=0; j<3; j++)
      mu_koll[k][j] = 0;
  }

  for(k=0; k<nz; k++) 
    for(j=0; j<3; j++)
      for(i=0; i<nz; i++)     
	mu_koll[k][j] += eigvec[i][k]*qy[i][j]*dipst[i];
    
  for(k=0; k<nz; k++){
    for(j=0; j<3; j++){
      cos_theta[k] += sym_ax[j]*mu_koll[k][j]; 
    }
    cos_theta[k] = SQ(cos_theta[k]);
  }
   
  for(k=0; k<nz; k++)
    ld_koeff[k] = 0.5*alpha[k]-1.5*cos_theta[k];

} 


//============================================================================


void sumKoeff(int i, int nz, double **eigvec, double **sumkoeff){

  // c_m^(M) = eigvec[m][M] ist m-te Komponente des M-ten Eigenvektors
  // m-tes Pigment, M-tes Exciton <=> c[Pigment][Exciton]=c[m][M]
  // kollektives Dipolmoment  mu_K = sum_i{ c_i^(K)*mu_i } mit
  // lokalen Dipolmomenten mu_i

  int j=0, k=0;

  for(k=0; k<nz; k++)
    for(j=0; j<nz; j++)
      sumkoeff[k][j] += SQ(eigvec[k][j]);

}

//============================================================================================

void spekDiff(int anz, double *spek_ab, double *spek_hb, double *spek_hbdif)
{
  // Berechnet homogenen HB-Differenzspektrum, delta alpha_h,k,(i)

  int j=0; 

  for(j = 0; j<anz; j++) //Berechnung des homogenen hb Differenzspektrums, delta alpha_h,k,(i)
  {
    spek_hbdif[j] = spek_hb[j] - spek_ab[j];
  }

}

//============================================================================================

void spekSum(int anz, double *spek, double *sumspek)
{
  // Aufsummieren der homogenen Spektren zu einem inhomog Spektrum

  int j=0; 

  for(j = 0; j<anz; j++) //Berechnung des homogenen hb Differenzspektrums, delta alpha_h,k,(i)
  {
    sumspek[j] += spek[j]/(1.0*N);
  }

}

//============================================================================================

void pigmExcDistrib(int nz, double *w_m0_sh, double **npig, double **nexc, double **nr_ev){

  //c[Pigment][Exciton]=c[m][M]

  int i=0, j=0, k=0;
  double x=0.0, g[nz][NEXC], gamha = 10.0;  //gamha=Gamma/2 = FHM
    
  for(k=0; k<nz; k++) for(j=0; j<NEXC; j++) g[k][j]=0.0;


  // Calculate exciton states pigment distribution
  // d_M(w)= < sum_M|c_m^(M)|^2 delta(w-w_M) >_dis

  for(i=0; i<nz; i++){

    for(k=0; k<nz; k++)
      for(j=0; j<NEXC; j++)
	g[k][j] = 0.0;

    for(k=0; k<nz; k++){
      for(j=0; j<NEXC; j++){
	x = SPEKMIN+(j+0.5)*EXCLENG;

	g[k][j] = SQ(nr_ev[i+1][k+1])*gamha/( SQ(x-w_m0_sh[k]) + SQ(gamha) );
      }
    }
    for(j=0; j<NEXC; j++){
      for(k=0; k<nz; k++){
	npig[i][j] += g[k][j];
      }
    }

  }//endfor_i


  // Calculate excitons:

  for(i=0; i<nz; i++){

    for(k=0; k<nz; k++)
      for(j=0; j<NEXC; j++)
	g[k][j] = 0.0;

    for(k=0; k<nz; k++){
      for(j=0; j<NEXC; j++){
	x = SPEKMIN+(j+0.5)*EXCLENG;

	g[k][j] = SQ(nr_ev[k+1][i+1])*gamha/( SQ(x-w_m0_sh[i]) + SQ(gamha) );
      }
    }
    for(j=0; j<NEXC; j++){
      for(k=0; k<nz; k++){
	nexc[i][j] += g[k][j];
      }
    }

  }//endfor_i


}


