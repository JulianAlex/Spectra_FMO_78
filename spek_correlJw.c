#include "spek_head.h"


void calcJw(double *jw){

/**************************************************************************
 void calcJw berechnet die Funktion J(w) 
 (entweder nach Renger & Marcus oder nach Pieper, siehe SPECDENS)

 Einheiten: J(w) hat die reziproke einheit von w, also wenn 
 [w] = fs-1  ==>  [J(w)] = fs
**************************************************************************/

  FILE *out_file; char *out_dat;

  double  w=0.0;
  int species=SPECIES, specdens=SPECDENS;

  int k=0, i=0;

  double fac = 2*5040;                   // =2*7! J(w)   [Renger & Marcus (2002), eq 2.17 ]
  double fac_pieper = 195.0119548*HRFAK; // Skalierungs-Parameter für Pieper J(w) so dass Huan-Rhys S = 1.3*HRFAK
  double fac_lognorm = 1.529412286*HRFAK; // Skalierungs-Parameter für LogNorm J(w) so dass Huan-Rhys S = 1.3*HRFAK

  double j_m, w_m[3], fw_gau[3], fw_lor[3], sig[3], s[3]; // For Pieper-J(w)

  out_dat="out_calcJw.dat";
  out_file = fopen(out_dat, "w");


  if(specdens == 1)       // Renger-Marcus spectral density
  {
    printf("Spectral Density from Renger and Marcus\n");  // SPECIES=%2d\n", species);

    for(k=0; k < NT/2; k++)
    {
      w = k*WSTEP;                // [fs-1]  

      jw[k] = CUB(w)*( S1/(fac*pow(W1,4))*exp(-sqrt(w/W1)) + S2/(fac*pow(W2,4))*exp(-sqrt(w/W2)));

      fprintf(out_file," %d  %lf  %lf  %lf\n", k, w, w/WZ2WFS, jw[k]); // J(w) stimmt mit Renger & Marcus ueberein. 
    }
  }


  else if(specdens == 2)  // Pieper spectral density //Pieper 2011, JPCB, Artcile I, S.4042
  { 
    if(species == 1)  // achtung, evt braucht man einen zweiten Skalierungs-Faktor!!! 
    {
      s[0]=0.45;  s[1]=0.20;  s[2]=0.16;              // a-WSCP
      w_m[0]=24.; w_m[1]=48.; w_m[2]=88.; 
      
      fw_gau[0]=25.; fw_gau[1]=20.; fw_gau[2]=24.;  // a-WSCP
      fw_lor[0]=25.; fw_lor[1]=50.; fw_lor[2]=60.;     

      printf("Spectral Density from Pieper, WSCP-A\n");
    }
    else if(species == 2)
    {
      s[0]=0.42;  s[1]=0.26;  s[2]=0.17;           // b-WSCP
      w_m[0]=23.; w_m[1]=46.; w_m[2]=82.; 
      
      fw_gau[0]=25.; fw_gau[1]=20.; fw_gau[2]=24.;  // b-WSCP
      fw_lor[0]=25.; fw_lor[1]=50.; fw_lor[2]=50.;   

      printf("Spectral Density from Pieper, WSCP-B\n");  
    }
    else
    {
      printf("WSCP A or B?! Check in spek_head.h!\n"); exit(8);
    }

    for(i=0; i<3; i++)
    {
      fw_gau[i] *= WZ2WFS;
      fw_lor[i] *= WZ2WFS;   
      w_m[i]   *= WZ2WFS;
      sig[i] = fw_gau[i]/FWHMSIG;
    }

    for(i=0; i<3; i++)
    {
      for(k=0; k<NT/2; k++) 
      {
	w = k*WSTEP;                // [fs-1]  
	    
	if(w <= w_m[i])
	{
	  j_m = fac_pieper*exp( -SQ(w-w_m[i])/(2*SQ(sig[i])) ); //Gauss-Kurve   
	}
	else if(w > w_m[i])
	{
	  j_m = fac_pieper*SQ(0.5*fw_lor[i])/( SQ(w-w_m[i]) + SQ(0.5*fw_lor[i]) ); //Lorentz-Kurve
	}
	else
	{     
	  printf("\nFehler in calcJw() !!\n"); exit(8);
	}

	jw[k] += s[i]*j_m;

	if(i==2)
	{	
	  fprintf(out_file," %d  %lf  %lf  %lf\n", k, w, w/WZ2WFS, jw[k]);
	} 
      }
    }


  }
  else if(specdens == 3)  // Log-Normal spectral density, Parameter: Kell 2013
  { 
    if(species == 1) // achtung, evt braucht man einen zweiten Skalierungs-Faktor!!! 
    {
      s[0]=0.; s[1]=0.; s[2]=0.;     // a-WSCP
      w_m[0]=0; w_m[1]=0; w_m[2]=0;       
      sig[0]=0; sig[1]=0; sig[2]=0;  // a-WSCP
           
      printf("Log-Normal Spectral Density, WSCP-A\n");
    }
    else if(species == 2)
    {
      s[0]=0.39; s[1]=0.23; s[2]=0.23;    // b-WSCP
      // w_m[0]=26; w_m[1]=51; w_m[2]=85;    // cm-1 // so kommen viel zu grosse Taus raus!!
      w_m[0]=52; w_m[1]=102; w_m[2]=170;    // cm-1 
      sig[0]=0.4; sig[1]=0.25; sig[2]=0.2;  // b-WSCP   

      printf("Log-Normal Spectral Density, WSCP-B\n");  
    }
    else
    {
      printf("WSCP A or B?! Check in spek_head.h!\n"); exit(8);
    }


    for(i=0; i<3; i++){
      w_m[i] *= WZ2WFS;          // cm-1 => fs-1

      jw[0] = 0;  

      for(k=1; k<NT/2; k++)
      {
	w = k*WSTEP;                // [fs-1]  
	
	j_m = fac_lognorm/(sqrt(TWOPI)*w)*s[i]/sig[i]*exp( -SQ(log(w/w_m[i]))/(2*SQ(sig[i])) );
	jw[k] += j_m;

	if(i==2)
	  fprintf(out_file," %d  %lf  %lf  %lf\n", k, w, w/WZ2WFS, jw[k]); 
      }
    }
  }

  else if(specdens == 4)       // modified Renger-Marcus SpecDens
  {
    printf("Modified spectral sensity from Renger and Marcus, Species=%2d\n", species);

    for(k=0; k < NT/2; k++)
    {
      w = k*WSTEP;                // [fs-1]  

      jw[k] = CUB(w)*( S31/(fac*pow(W31,4))*exp(-sqrt(w/W31)) + S32/(fac*pow(W32,4))*exp(-sqrt(w/W32)) + S33/(fac*pow(W33,4))*exp(-sqrt(w/W33)) );

      fprintf(out_file," %d  %lf  %lf  %lf\n", k, w, w/WZ2WFS, jw[k]); 
    }
  }

  else
  {
    printf("Spectral Density 1, 2, 3 or 4 ?! Check in spek_head.h!\n"); exit(8);
  }


  fclose(out_file);

}


//====================================================================================================================


void calcGt(fftw_complex *gt, double *jw, fftw_complex *in, fftw_complex *out, fftw_plan plan){

/**************************************************************************
 void calcGt berechnet die funktion G(t) = G(t,T) mit Hilfe einer
 Fouriertransformation (FFT mit num.recipes routine four1) von w nach t
 (ruecktrans, -1) = Int{exp(iwt)dw}

 G(t) = Int_0^Infty dw {(1+n(w))J(w)e^-iwt + n(w)J(w)e^iwt}
      = FFT-1{ f(w) } 

 mit f(w)=(1+n(w))*J(w) fuer w>0 und f(w)=n(|w|)*J(|w|) fuer w<0 

 w = 2*PI*f = 2*PI/t
**************************************************************************/

  FILE *out_file; char *out_dat;

  double  w=0.0, nj_pos=0.0, nj_neg=0.0, nvib=0.0;
  int k=0;

  for(k=0; k<NT; k++)
  {
    in[k] = 0.0;
    out[k] = 0.0;
  }

  for(k=1; k<NT/2; k++)         // j[0]=0;  nvib(w=0) => nan!!
  {
    w = k*WSTEP;                // [fs-1]  

    nvib = 1.0/( exp(w*HBDKB/TEMP)-1.0 );    // nvib(w=0) => nan!!
    
    nj_pos = jw[k]*(1.0+nvib);
    nj_neg = jw[k]*nvib;        // [fs]
    
    in[k]        = nj_pos;  // Real-Teil positive Frequenzen
    in[NT-(k+1)] = nj_neg;  // Real-Teil negative Frequenzen
  }


  fftw_execute(plan);   // FFT  h_k = 1/N*sum{H_n e^-ikn/N}


  out_dat="out_gt.dat";
  out_file = fopen(out_dat, "w");
 
  for(k=0; k<NT; k++)
  {
    gt[k] = TWOPI/(DELTA*NT)*out[k];   // gt complex

    fprintf( out_file," %d  %lf  %lf  %lf\n", k, k*DELTA, creal(gt[k]), cimag(gt[k]) ); // t=k*Delta [fs]
  }

  fclose(out_file);


}

//==========================================================================================================


void calcCt(fftw_complex *ct, double *jw, fftw_complex *in, fftw_complex *out, fftw_plan plan){

/**************************************************************************
 void calcCt berechnet die Funktion C(t) mit Hilfe einer
 Fouriertransformation (FFT mit num.recipes routine four1) von w nach t
 (ruecktrans, -1) = Int{exp(iwt)dw}

 C(t) = Int_0^Infty dw w^2 {(1+n(w))J(w)e^-iwt + n(w)J(w)e^iwt}
      = FFT-1{ w^2*f(w) }  

 Achtung: ct[] = gt[0...2NT-1]  

 mit f(w)=(1+n(w))*J(w) fuer w>0 und f(w)=n(|w|)*J(|w|) fuer w<0 

 w = 2*PI*f = 2*PI/t 

 Fuer 77K gleiches C(t) wie Reng&Marc 2002, Abb 3
**************************************************************************/

  FILE *out_file; char *out_dat;
  out_dat="out_calcCtJw.dat";
  out_file = fopen(out_dat, "w");

  double  w=0.0, nj_pos=0.0, nj_neg=0.0, nvib=0.0;
  int k=0;

  for(k=0; k<NT; k++)
  {
    in[k] = 0.0;
    out[k] = 0.0;
  }

 
  for(k=1; k<NT/2; k++)  // j[0]=0,  nvib(w=0) => nan!! 
  {
    w = k*WSTEP;   // [fs-1]  
    
    nvib = 1.0/( exp(w*HBDKB/TEMP)-1.0 );   // nvib(w=0) => nan!! 

    nj_pos = jw[k]*SQ(w)*(1.0+nvib);  // [fs-1]
    nj_neg = jw[k]*SQ(w)*nvib;      

    in[k]        = nj_pos;  // Real-Teil positive Frequenzen
    in[NT-(k+1)] = nj_neg;  // Real-Teil negative Frequenzen

    fprintf(out_file," %lf  %lf  %lf  %lf  %lf\n", w, w/WZ2WFS, jw[k], nj_pos, nj_neg ); // J(w) stimmt mit Renger & Marcus ueberein. 
  }
  fclose(out_file);

  

  fftw_execute(plan);   // inverse FFT,  FT{in} = C(t) ------------------------


  out_dat="out_ct.dat";
  out_file = fopen(out_dat, "w");
 
  for(k=0; k<NT; k++)
  {
    ct[k] = TWOPI/(DELTA*NT)*out[k];   // [ct]=fs-2

    fprintf(out_file," %d  %lf  %g  %g\n", 
    	    k, k*DELTA, creal(ct[k])*SQ(HBAR*1000), cimag(ct[k])*SQ(HBAR*1000) ); 
  }

  fclose(out_file);  

}  // C(t) stimmt mit Renger & Marcus ueberein (fuer T=77K !!).


// ===============================================================================================


void calcCw(fftw_complex *ct, fftw_complex *cw, double *cw_re, double *cw_im, fftw_complex *in, fftw_complex *out, fftw_plan plan){

/**************************************************************************
 void calcCW berechnet die funktion ~C(w) als halbseitige
 Fouriertransformation von C(t)
 H(w)=Int{C(t)e^iwt}dt, Hintrans +1) 
**************************************************************************/
  // DELTA  [fs]  
  // WSTEP  [fs-1] 

  FILE *out_file; char *out_dat;
 
  int k=0;

  for(k=0; k<NT; k++)
  {
    in[k] = 0.0;
    out[k] = 0.0;
  }


  for(k = 0; k < NT/2; k++)   // halbseitige FT, dh Eintraege fuer k > NT/2 sind null!
  {
    in[k] = ct[k];   // [fs-2]  
  }
 

  fftw_execute(plan);   // forward-FFT,  C(w)=Int{ C(t) exp(-iwt) } dt 

  
  for(k=0; k<NT; k++)   // wegen omegaShift wird cw in cm-1 benoetigt
  {
    cw[k] = DELTA/WZ2WFS*out[k];  
  }
  // [out]=fs-2, [DELTA]=fs => [nr_data*DELTA]=fs-1  
  // [out*DELTA]/WZ2FS = [cw] = cm-1


  for(k=NT/2; k<NT; k++)
  {
    cw_re[k-NT/2] = creal(cw[k]);    //cw_im[k] fuer k = 0...NT/2 negativ, NT/2...NT positiv 
    cw_im[k-NT/2] = cimag(cw[k]);    
  }
  for(k=0; k<NT/2; k++)
  {
    cw_re[k+NT/2] = creal(cw[k]); 
    cw_im[k+NT/2] = cimag(cw[k]);      // [~C(w)] = cm-1
  }

 
  out_dat="out_cw.dat";
  out_file = fopen(out_dat, "w");

  for(k=0; k<NT; k++)
    fprintf(out_file," %d  %lf  %lf  %lf\n", k, (k-NT/2)*WSTEP, cw_re[k], cw_im[k]); 

  fclose(out_file);
 
}
 
