#include "spek_head.h"
//
// Parametrisierte Version
//
// Julian Adolphs
//
// Funktion "lineshape" berechnet die Linienform nach dynamischer Theory
//
// Die Fouriertransformation in lineshape ist eine von t nach w (hintrans, 1)
// Fouriertransformationen mit FFTW
// Einheiten:  [zeroline] = fs
//
// Das gamma-tau-Array: 
//
// gamma = [DGAM_0, DGAM_0+DGAM, DGAM_0+2*DGAM,..., DGAM_0+(NGAM-1)*DGAM] 
// 
// tau   = [DTAU_0, DTAU_0+DTAU,..., DGAM_0+(NGAM-2)*DGAM, TAU_0]  
//
// dh. gamma lauft einfach in DGAM-Schritten von gamma_min bis gamma_max, 
// tau laeuft ebenfalls von tau_min bis tau_max in DTAU-Schritten, ABER 
// der letzte Eintrag springt dann auf den Wert TAU_0. da zwischen tau_max und TAU_0
// eine riesige Luecke ist und es keinen Sinn macht, das Array mit nicht benoetigten 
// Werten vollzuschreiben!! 

void lineshapeArray(fftw_complex *gt, float ***sideband, float ***zeroline, 
		    fftw_complex *in, fftw_complex *out, fftw_plan plan){

  //  FILE *out_file; char *out_dat;

  int k=0, m=0, n=0;
  double t=0.0, re=0.0, im=0.0;

  //double wstep = WSTEP, ww=0.0;  // Schrittweite [fs-1] ??

  double gam[NGAM], tau[NTAU];   // gamma-tau-Array

  for(m=0; m<NGAM; m++)
    gam[m] = m*DGAM+DGAM_0;

  for(m=0; m<NTAU-1; m++)
  {
    tau[m] = m*DTAU+DTAU_0;    //in fs
  }
  tau[NTAU-1] = TAU_0; // letzter array-eintrag enhält maximales tau=TAU_0


  double reFT[NT];

  for(n=0; n<NT; n++) 
    reFT[n]=0.0;


  
  for(m=0; m<NGAM; m++)   // ----- parametrisier-loop ----------
  {
    for(n=0; n<NTAU; n++)
    {

      // zero vibrational quanta line (including neg freq) with FFT
      
      for(k=0; k<NT/2; k++)
      { 
	t = DELTA*k; 
	in[k] = exp( -t/tau[n] );
      }   
      for(k=NT/2; k<NT; k++)
      { 
	t = DELTA*(NT-k); 
	in[k] = exp( -t/tau[n] );
      }   
          

      fftw_execute(plan);      // FFT Fourier-Transform  (vorwaerts) 


      for(k=0; k<NT; k++)            
      	reFT[k] = TWOPI*DELTA*creal(out[k]); // Realteil der FT


      for(k=0; k<NT/2; k++)
	zeroline[m][n][k] = reFT[k+NT/2];  

      for(k=NT/2; k<NT; k++)
	zeroline[m][n][k] = reFT[k-NT/2];     

    

      // vibrational sideband via Fouriertransform: 
    
      for(k=0; k<NT/2; k++)
      {
	t = DELTA*k;
	re = gam[m]*creal(gt[k]);    //G_M(t)=gam_MM*G(t)           
	im = gam[m]*cimag(gt[k]);
	in[k] = exp(-t/tau[n])*( (exp(re)*cos(im)-1) + I*exp(re)*sin(im) );     
      }

      for(k=NT/2; k<NT; k++)
      {
	t = DELTA*(NT-k);
	re = gam[m]*creal(gt[k]);
	im = gam[m]*cimag(gt[k]);
	in[k] = exp(-t/tau[n])*( (exp(re)*cos(im)-1) + I*exp(re)*sin(im) );
      }
   
 
      fftw_execute(plan);   // --- FFT (vorwaerts) -------------------


      for(k=0; k<NT; k++)           
	reFT[k] = TWOPI*DELTA*creal(out[k]);   // Real part of FT  
      

      for(k=0; k<NT/2; k++)
	sideband[m][n][k] = reFT[k+NT/2];  
    
      for(k=NT/2; k<NT; k++)
	sideband[m][n][k] = reFT[k-NT/2];     
       
    } // end n-loop

  } // end m-loop

}

 

/*************************************************************************************

Tests: Verhältnis (bei T=0K) SB/ZL = 1.159 = exp(G_M(0))-1 = 1.159   ok!  

       Integral uber das Spektrum (fast) nicht temperaturabhaengig (< 1% Abweichung) 

     
// es gilt zeroline = 2*markov, das laesst sich ueberpruefen, weil fuer einen
// lokalen uebergang gilt: Int(SB)/Int(ZL)=exp(gamma_MM*S)-1 (siehe unten)




   if( n==1 )
      {
	out_dat="out_lineshape.dat";
	out_file = fopen(out_dat, "w");
	for(k=0; k<NT; k++){ 
	  ww = wstep*k;
	  fprintf(out_file," %d  %lf  %lf  %lf\n", k, ww,
		  sideband[n][k], zeroline[n][k]); 
	}
	fclose(out_file); 
      }
 

    if( n==0 )
      {
	out_dat="out_zeroline.dat";
	out_file = fopen(out_dat, "w");
	for(k=0; k<NT; k++){ 
	  fprintf(out_file," %d  %lf\n", k, zeroline[n][k] ); 
	}
	fclose(out_file); 
      }



    if( n==0 )
      {
	out_dat="out_sideband.dat";
	out_file = fopen(out_dat, "w");
	for(k=0; k<NT; k++){ 
	  fprintf(out_file," %d  %lf  %lf\n", k, sideband[n][k]); 
	}
	fclose(out_file); 
      }



// for(n=0; n<nz; n++) 
//   printf("G_M(0) = %lf, G(t=0)= %lf    Remark: G(t=0,T=0) = S  \n", gam[n][n]*creal(gt[0]), creal(gt[0]) ); printf("\n");



************/
