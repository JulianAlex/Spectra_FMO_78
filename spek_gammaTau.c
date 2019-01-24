#include "spek_head.h"


// Re{G_M(0)} = gmt[m][NT]  fuer T=0K => Re{G_M(0)}=S  Huang Rhys
// Im{G_M(0)} = 0 !!!

// Spectral density J(w): jw[j]
// Number of vibrational quanta n(w): nvib[j]
// n(w)=1/(exp(hw/kT)-1),  hw/kT --> facn*w(1/cm)/T

//===========================================================================

void gammaTau(int nz, double **rij, double *eigval, double **eigvec, 
	      double **gamma, double *tau, double *cw_re, double *jw){

  // FILE *out_file; char *out_dat;

  int j=0, k=0, m=0, mm=0, jj=0, kk=0;

  double cw_re_wmk[nz][nz], w=0.0, tau_rec[nz];

  // WSTEP [fs-1]

  for(j=0; j<nz; j++)
  {
    for(k=0; k<nz; k++)
    {
      gamma[j][k]     = 0.0; 
      cw_re_wmk[j][k] = 0.0;
    }
    tau[j]     = 0.0;
    tau_rec[j] = 0.0;
  }


  // cw_re_wmk lieber direkt statt aus der FFT berechnen! Ist genauer. 
  // Gibt doch eine relative grosse Abweichung fuer Tau, wenn man cw_re[m] verwendet...


  // Berechne ~C_Re(w) = pi w^2 { (1+n(w))J(w) + n(-w)J(-w) }   R&M (4.15)

  double nwjw=0.0, nvib=0.0;   

  for(m=0; m<nz; m++) 
  {
    for(k=0; k<nz; k++)
    {
      w = (eigval[m] - eigval[k])*WZ2WFS;   // [fs-1]  

      if(w > 0)
      {
	nvib = 1.0/( exp( w*HBDKB/TEMP )-1.0 ); 
	jj   = round(w/WSTEP); 
	nwjw = jw[jj]*SQ(w)*(1 + nvib); // [fs-1] 
      }	
      if(w < 0)
	{
	  nvib = 1.0/( exp( -w*HBDKB/TEMP )-1.0 ); 
	  jj   = round(ABS(w)/WSTEP); 
	  nwjw = jw[jj]*SQ(w)*nvib; 
	}
      if(w == 0)
      {
	jj   = 0;
	nwjw = 0;
      }
      cw_re_wmk[m][k] = PI*nwjw; 
      
      // printf("wmk, C^Re(w)  %10.1lf  %e  %e\n", w/WZ2WFS, cw_re_wmk[mm][kk], jott);
    }
  }


  for(mm=0; mm<nz; mm++)         // eq 4.9  &  eq 4.11
    for(kk=0; kk<nz; kk++)       // mm und kk = exciton state
      for(m=0; m<nz; m++)        // m = Pigment
	gamma[mm][kk] += SQ(eigvec[m][mm])*SQ(eigvec[m][kk]); 


  /* wenn Korrelations(-radius) erwuenscht, Aenderungen in gammaTau.c hier einkommentieren!!!     

  for(mm=0; mm<nz; mm++)         // eq 4.9  &  eq 4.11
    for(kk=0; kk<nz; kk++)       // mm und kk = exciton state
      for(m=0; m<nz; m++)        // m = Pigment
	for(n=0; n<nz; n++)
	  gamma[mm][kk] += exp(-rij[m][n]/RCORR)*eigvec[m][mm]*eigvec[m][kk]*eigvec[n][kk]*eigvec[n][mm]; 
          
   // Einheiten: [eigvec] = [c_m^(M)] = 1,  [gamma] = 1     
  */



  for(mm=0; mm<nz; mm++)    // eq 4.22  &  eq 4.12
    for(kk=0; kk<nz; kk++)
      tau_rec[mm] += gamma[mm][kk]*cw_re_wmk[mm][kk];  // tau_rec = 1/tau,  [tau_rec]=fs-1


  for(mm=0; mm<nz; mm++)
  {
    tau[mm] = 1.0/(tau_rec[mm]*LIFETIMEBROADENING + 1.0/TAU_0);    // [fs]
  }

  /*
  out_dat="out_gammaTau.dat";
  out_file = fopen(out_dat, "a");
  //for(kk=0; kk < nz; kk++)
  //  fprintf(out_file, "gammTau.c  %3d   gamma = %lf,  tau= %8.1lf fs,   %6.2lf cm-1\n", 
  //	   kk, gamma[kk][kk], tau[kk], 1./tau[kk]/WZ2WFS);
  //  fprintf(out_file, "gammTau.c  %3d   gamma = %lf,  tau= %8.1lf fs,   %6.2lf cm-1\n", 
  //	   1, gamma[1][1], tau[1], 1./tau[1]/WZ2WFS);
  fprintf(out_file, " %8.1lf \n", tau[1]);  // tau [fs]
  fclose(out_file);
  */

}

// ======================================================================================================0

