/****************************************************************************

Julian Adolphs 

Programm zu Berechnung der inhomogenen linearen Spektren (Absorption, CD, LG) 
von FMO (mit 7 und 8 Pigmenten in vorgegebenen Mischungsverhaeltnis) 
in dynamischer Theorie mit einer MonteCarlo-Methode.

Input: X-Ray-Koordinaten, Site-Energies (spek.h)

Output: out_spec_OD_nonMark.dat, out_spec_CD_nonMark.dat, out_spec_LD_nonMark.dat

Files die zum Programm "spektrum" gehoeren:  Makefile, spek_*.c, spek_*.h

Fourier-Transformationen mit FFTW-Routinen 

Diagonalisierung der Hamilton-Matrix mit "Eispack":  eispack.c, eispack.h

Zufallszahlen:  normal.h, normal.c 

Theorie: T. Renger & R.A. Marcus (2002), J. Chem. Phys 116 (22), 9997.  
On the relation of protein dynamics and exciton relaxation in pigment-protein complexes
File: Reng02b.pdf

****************************************************************************/
  

#include "spek_head.h"   // Makros & Definitionen


int main(int argc, char *argv[ ]) 
{

  clock_t prgstart, prgende;  // Laufzeitmessung
  prgstart = clock();         // CPU-Zeit zu Beginn des Programmes

  int i=0, k=0, n=0, nz=0, nv=0, trim_nz=0, hn=0, nrow=0, ncol=0; // pigm_no=0; 

  int bclnum = 8;        // 7 or 8 pigment variant of FMO

  int nmc_8 = N*BCL_8;   // Anteil MC-steps mit 8 pigmenten

  int h_var = 17;        // seed for random number generator
  int *seed  = &h_var;

  int couplings = COUPLINGS;   // Kopplungen berechnen oder einlesen?!

  //  int species   = SPECIES;     // Chl_ or Chl_b-WSCP or Tep-FMO or Aest-FMO  
  //  int dyntheo   = DYNTHEO;     // Dynamische Theorie mit Markov oder nonMarkov

  char *in_dat = "infile_koord.dat";
  n = countInput(in_dat);

  nz = n/NNATM;                  // anzahl der chlorophylle

  trim_nz = 3*nz;
  hn = SQ(nz);
  printf("# Anzahl N-Atome:  %d ,  Anzahl Chl:  %d,  Temperatur:  %4.1lf K\n", n, nz, TEMP);
  
  double wstep = WSTEP/WZ2WFS;  //Schrittweite w-Achse [cm-1]  
  int    anz   = (int) (SPEKDIFF/wstep + 1);    
  printf("# wstep:  %lf cm-1, anz:  %d\n", wstep, anz);


 // ------------- Declare & Allocate vectors  ---------------------------------

  int   *num      = calloc( n, sizeof(int));
  int   *trim_num = calloc( trim_nz, sizeof(int) );


  double *x = calloc( n, sizeof(double));
  double *y = calloc( n, sizeof(double));
  double *z = calloc( n, sizeof(double));

  double *sym_ax =  calloc( 3, sizeof(double));    // Symmetrie-Achse fuer LD-Spek

  double *dipst   = calloc( nz, sizeof(double));
  double *site    = calloc( nz, sizeof(double));   // Site Energie Mittelwerte 
  double *sigma   = calloc( nz, sizeof(double));
  double *alpha   = calloc( nz, sizeof(double));   // alpha=|mu|^2 dipst^2 [in Debye^2] 
  double *rotst   = calloc( nz, sizeof(double));   // rotational strength for CD 
  double *ldcff   = calloc( nz, sizeof(double));   // LD-coefficient for LD 
  double *tau     = calloc( nz, sizeof(double));
  double *w_m0_sh = calloc( nz, sizeof(double));

  double *jw    = calloc( NT/2, sizeof(double));   // J(w)   jw[] = jw[0] ... jw[NT/2-1]

  double *cw_re = calloc( NT, sizeof(double));
  double *cw_im = calloc( NT, sizeof(double));

  double *spek_ab       = calloc( anz, sizeof(double)); // Homogenes Absorptions-Spektrum   
  double *spek_sb       = calloc( anz, sizeof(double)); // Homogenes Absorptions-Spektrum nur der SeitenBande
  double *spek_zl       = calloc( anz, sizeof(double)); // Homogenes Absorptions-Spektrum nur der ZeroVibLine
  double *spek_cd       = calloc( anz, sizeof(double)); // Homogenes CD-Spektrum 
  double *spek_ld       = calloc( anz, sizeof(double)); // Homogenes LD-Spektrum 
  double *spek_cdmk     = calloc( anz, sizeof(double)); // Homogenes CD-Spektrum in Markov-Approx
  double *spek_ldmk     = calloc( anz, sizeof(double)); // Homogenes LD-Spektrum in Markov-Approx
  double *sumspek_ab    = calloc( anz, sizeof(double)); // Inhomogenes Absorptionsspektrum (=aufsummiert)   
  double *sumspek_sb    = calloc( anz, sizeof(double)); // Inhomogenes AbsSpek der SB (=aufsummiert)
  double *sumspek_zl    = calloc( anz, sizeof(double)); // Inhomogenes AbsSpek der ZL (= Markov-Approx)
  double *sumspek_cd    = calloc( anz, sizeof(double)); // Inhomogenes CD-spektrum (=aufsummiert)
  double *sumspek_ld    = calloc( anz, sizeof(double)); // Inhomogenes LD-spektrum (=aufsummiert)
  double *sumspek_cdmk  = calloc( anz, sizeof(double)); // Inhomogenes CD-spektrum in Markov-Approx
  double *sumspek_ldmk  = calloc( anz, sizeof(double)); // Inhomogenes LD-spektrum in Markov-Approx

  // Achtung: fftw_malloc initialisiert den Speicher nicht! 
  fftw_complex *gt; gt = (fftw_complex*) fftw_malloc( NT * sizeof(fftw_complex) );   // G(t)   
  fftw_complex *ct; ct = (fftw_complex*) fftw_malloc( NT * sizeof(fftw_complex) );   // C(t)   
  fftw_complex *cw; cw = (fftw_complex*) fftw_malloc( NT * sizeof(fftw_complex) );   //~C(w) 

  for(i=0; i<NT; i++)
  {
    gt[i] = 0.0;
    ct[i] = 0.0;
    cw[i] = 0.0;
  }

  fftw_complex *in, *out;
  fftw_plan plan_bwd, plan_fwd;

  in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NT); //allocating memory
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NT);

  for(k=0; k<NT; k++)
  {
    in[k]  = 0.0;
    out[k] = 0.0;
  }

  //  plan = fftw_plan_dft_1d(NT, in, out, +1, FFTW_MEASURE);  // -1/+1 direction of FT, FFTW_MEASURE accuracy
  plan_bwd = fftw_plan_dft_1d(NT, in, out, +1, FFTW_MEASURE);  // backward-transform
  plan_fwd = fftw_plan_dft_1d(NT, in, out, -1, FFTW_MEASURE);  // forward-transform


  //---- Matrizen -----------------------------------------------------------------------------------
  // Im FFTW-Manual wird von dieser Art Arrays zu konstruieren abgeraten, Schmaranz empfiehlt es so...

  double **mg = calloc( nz, sizeof(double*) );
  double **na = calloc( nz, sizeof(double*) );
  double **nb = calloc( nz, sizeof(double*) );
  double **nc = calloc( nz, sizeof(double*) );
  double **nd = calloc( nz, sizeof(double*) );
  double **qx  = calloc( nz, sizeof(double*) );
  double **qy   = calloc( nz, sizeof(double*) );
  double **vab   = calloc( nz, sizeof(double*) );
  double **rij   = calloc( nz, sizeof(double*) );
  double **gam  = calloc( nz, sizeof(double*) );
  double **j_w = calloc( nz, sizeof(double*) );     // Spectral Density J(w)
  double **n_w = calloc( nz, sizeof(double*) );     // No of vib.quanta n(w)

  nrow = nz;
  while(nrow--)
  {
    mg[nrow] = calloc(3, sizeof(double));
    na[nrow] = calloc(3, sizeof(double));
    nb[nrow] = calloc(3, sizeof(double));
    nc[nrow] = calloc(3, sizeof(double));
    nd[nrow] = calloc(3, sizeof(double));
    qx[nrow] = calloc(3, sizeof(double));
    qy[nrow] = calloc(3, sizeof(double));
    vab[nrow] = calloc(nz, sizeof(double));
    rij[nrow] = calloc(nz, sizeof(double));
    gam[nrow] = calloc(nz, sizeof(double));
    j_w[nrow] = calloc(nz, sizeof(double));
    n_w[nrow] = calloc(nz, sizeof(double));
  }

  double **trim_r = calloc( trim_nz, sizeof(double*));
  nrow = trim_nz;
  while(nrow--)
    trim_r[nrow] = calloc(3, sizeof(double));

  double **npig = calloc( nz, sizeof(double*) );
  double **nexc = calloc( nz, sizeof(double*) );
  nrow = nz;
  while(nrow--)
  {
    npig[nrow]  = calloc( NEXC, sizeof(double));
    nexc[nrow]  = calloc( NEXC, sizeof(double));
  }

  double **delta = calloc( hn, sizeof(double*));
  nrow = hn;
  while(nrow--)
    delta[nrow] = calloc(3, sizeof(double));
  
  double **arrayspek_sb = calloc( nz, sizeof(double*) );
  double **arrayspek_zl = calloc( nz, sizeof(double*) );

  nrow = nz;
  while(nrow--)
  {
    arrayspek_sb[nrow] = calloc( 2*anz, sizeof(double) ); //calloc( 4*anz+9, sizeof(double) ); 
    arrayspek_zl[nrow] = calloc( 2*anz, sizeof(double) );
  }

 
  // float statt double halbiert den benötigten Arbeitsspeicher
  float ***sideband = calloc( NGAM, sizeof(float**) );          //[gamma][tau][w]
  float ***zeroline = calloc( NGAM, sizeof(float**) );          //[gamma][tau][w]
  nrow = NGAM;
  while(nrow--){
    sideband[nrow] = calloc( NTAU, sizeof(float*) );
    zeroline[nrow] = calloc( NTAU, sizeof(float*) );
    ncol = NTAU;
    while(ncol--){
      sideband[nrow][ncol] = calloc( NT, sizeof(float));
      zeroline[nrow][ncol] = calloc( NT, sizeof(float)); 
    }
  }

  
  double *hamilt    = calloc( nz*nz, sizeof(double));    // Hamilton-Matrix
  double *ev_h      = calloc( nz*nz, sizeof(double));    // Eigenvektoren (1D) output of the eispack-routines
  double *eigval    = calloc( nz, sizeof(double));       // Eigenwerte    

  double **eigvec    = calloc(nz, sizeof(double*));    // Eigenvektoren (2D)
  nrow = nz;
  while(nrow--)
  {
    eigvec[nrow]    = calloc(nz, sizeof(double));
  }


 //======= Begin Hauptprogramm "spektrum" ====================================

  in_dat="infile_koord.dat";        
  readCoordinates(num, x, y, z, n, in_dat);  
  parameter(site, sigma, dipst, bclnum, nz);        
  koord(nz, x, y, z, na, nb, nc, nd);
  pigmCenters(nz, mg, nb, nd, delta, rij);

  in_dat="infile_koord_trim.dat";   
  readTRImerCoord(trim_num, trim_r, trim_nz, in_dat);
  trimSymAx(trim_nz/3, trim_r, sym_ax); 


  if(couplings == 1) 
  {
    dipolmom(nz, na, nb, nc, nd, qx, qy);
    interaction(nz, qy, delta, vab, rij, dipst);
    printf("PointDipoleApproximation\n"); printf("\n");
  }
  else if(couplings == 2) 
  {  
    in_dat="infile_couplings.dat";
    nv = countInput(in_dat);
    dipolmom(nz, na, nb, nc, nd, qx, qy);
    readCouplings(nv, nz, vab, in_dat);
    printf("Couplings from File\n"); printf("\n");
  }
  else
  {
    printf("Exit with couplings = %d", couplings);
    exit(8);
  }
  if ( (E_BURN <= SPEKMIN) || (E_BURN >= SPEKMAX ) )
  {
    printf("\nFehler: Brennfrequenz außerhalb SPEKMIN und SPEKMAX\n");
    exit(8);
  }

//---------------------------------------------------------------------------------------------------------------

  calcJw(jw);                                       // Berechnung von J(w) eq 2.17,  [J(w)] = fs
  calcGt(gt, jw, in, out, plan_fwd);                // Berechnung von G(t) mit FFT,  [G(t)] = 
  calcCt(ct, jw, in, out, plan_fwd);                // Berechnung von C(t) mit FFT,  [C(t)] = fs-1
  calcCw(ct, cw, cw_re, cw_im, in, out, plan_bwd);  // Berechnung von C(w) mit FFT,  [C(w)] = cm-1

  
  lineshapeArray(gt, sideband, zeroline, in, out, plan_bwd);  // FastFourierTransf


  for(i=0; i<nmc_8; i++)  // BEGIN  Monte-Carlo-Loop 
  {
    calcLinearSpectrum( nz, i, anz, site, sigma, seed, vab, hamilt, eigval, eigvec, ev_h, qy, alpha, rotst, ldcff, 
			delta, dipst, rij, gam, tau, cw_re, cw_im, jw, gt, w_m0_sh, sideband, zeroline, 
			spek_ab, sumspek_ab, spek_sb, arrayspek_sb, sumspek_sb, spek_zl, arrayspek_zl, sumspek_zl, 
			spek_cd, sumspek_cd, spek_ld, sumspek_ld, spek_cdmk, sumspek_cdmk, spek_ldmk, sumspek_ldmk, 
			sym_ax, in, out, plan_bwd );   
  } // END 8 
    
  bclnum = 7;  
  parameter(site, sigma, dipst, bclnum, nz);
  dipolmom(nz, na, nb, nc, nd, qx, qy);
  interaction(nz, qy, delta, vab, rij, dipst);

  for(i=nmc_8; i<N; i++)  // BEGIN  7 Pigments
  {
    calcLinearSpectrum( nz, i, anz, site, sigma, seed, vab, hamilt, eigval, eigvec, ev_h, qy, alpha, rotst, ldcff, 
			delta, dipst, rij, gam, tau, cw_re, cw_im, jw, gt, w_m0_sh, sideband, zeroline, 
			spek_ab, sumspek_ab, spek_sb, arrayspek_sb, sumspek_sb, spek_zl, arrayspek_zl, sumspek_zl, 
			spek_cd, sumspek_cd, spek_ld, sumspek_ld, spek_cdmk, sumspek_cdmk, spek_ldmk, sumspek_ldmk, 
			sym_ax, in, out, plan_bwd );  
  } //END 7 Pigmets
    
  outputSpectrum(sumspek_ab, anz,   "out_spec_OD_nonMark.dat");   // absorption dynTheory nonMarkov
  outputSpectrum(sumspek_cd, anz,   "out_spec_CD_nonMark.dat");   // circular-dichrois dynTheory nonMarkov
  outputSpectrum(sumspek_ld, anz,   "out_spec_LD_nonMark.dat");   // linear-dichrois nonMarkov
	
  prgende=clock(); //CPU-Zeit am Ende des Programmes
  printf("\nLaufzeit %.2f Sekunden\n\n",(float)(prgende-prgstart) / CLOCKS_PER_SEC); 

// ===== END Hauptprogramm =========================================================================


// --- free allocated memory:

  free(num);  num  = NULL;

  free(x); x = NULL;
  free(y); y = NULL;
  free(z); z = NULL;

  free(dipst); dipst = NULL;
  free(site);  site = NULL;
  free(sigma); sigma = NULL;
  free(alpha); alpha = NULL;
  free(rotst); rotst = NULL;  
  free(ldcff); ldcff = NULL;
  free(tau);   tau = NULL;
  free(w_m0_sh); w_m0_sh = NULL;

  free(jw); jw  = NULL;
 
  free(cw_im); cw_im = NULL;
  free(cw_re); cw_re = NULL;


  free(spek_ab);    spek_ab = NULL; 
  free(spek_sb);    spek_sb = NULL; 
  free(spek_zl);    spek_zl = NULL; 
  free(spek_cd);    spek_cd = NULL; 
  free(sumspek_ab); sumspek_ab = NULL;
  free(sumspek_sb); sumspek_sb = NULL;
  free(sumspek_zl); sumspek_zl = NULL;
  free(sumspek_cd); sumspek_cd = NULL;


  // FFTW-complex-numbers

  fftw_free(ct); ct = NULL; 
  fftw_free(gt); gt = NULL;
  fftw_free(cw); cw = NULL;


// free vectors and matrices -------------------------------------------------------------


  nrow = nz;
  while(nrow--)
  {
    free(mg[nrow]);
    free(na[nrow]);
    free(nb[nrow]);
    free(nc[nrow]);
    free(nd[nrow]);
    free(qx[nrow]);
    free(qy[nrow]);
    free(vab[nrow]);
    free(rij[nrow]);
    free(gam[nrow]);
    free(j_w[nrow]);
    free(n_w[nrow]);
    free(npig[nrow]);
    free(nexc[nrow]); 
    free(arrayspek_sb[nrow]); 
    free(arrayspek_zl[nrow]);
  }
  free(mg); mg = NULL;
  free(na); na = NULL;
  free(nb); nb = NULL;
  free(nc); nc = NULL;
  free(nd); nd = NULL;
  free(qx); qx = NULL;
  free(qy); qy = NULL;
  free(vab); vab = NULL;
  free(rij); rij = NULL;
  free(gam); gam = NULL;
  free(j_w); j_w = NULL;
  free(n_w); n_w = NULL;
  free(npig); npig = NULL;
  free(nexc); nexc = NULL;
  free(arrayspek_sb); arrayspek_sb = NULL;
  free(arrayspek_zl); arrayspek_zl = NULL;


  nrow = hn;
  while(nrow--)
  {
    free(delta[nrow]);
  }
  free(delta); delta = NULL;


  nrow = NGAM;
  while(nrow--)
  {
    ncol = NTAU;
    while(ncol--)
    {
      free(sideband[nrow][ncol]);
      free(zeroline[nrow][ncol]); 
    }
    free(sideband[nrow]);
    free(zeroline[nrow]);   
  }
  free(sideband);
  free(zeroline);   
  

  free(eigval); eigval = NULL;
  free(ev_h);   ev_h   = NULL;
  free(hamilt); hamilt = NULL;

  nrow = nz;
  while(nrow--)
  {
    free(eigvec[nrow]);
  }
  free(eigvec); eigvec = NULL;

  

  // FFTW-plans
  
  fftw_destroy_plan(plan_bwd);
  fftw_destroy_plan(plan_fwd); 
  fftw_free(in); 
  fftw_free(out);
  

  return (0);

}

