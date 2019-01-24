#include "spek_head.h"

void calcLinearSpectrum(int nz, int i, int anz, double *site, double *sigma, int *seed, double **vab, double *hamilt, 
			double *eigval, double **eigvec, double *ev_h, double **qy, double *alpha, double *rotst, double *ldcff, 
			double **delta, double *dipst, double **rij, double **gam, double *tau, 
			double *cw_re, double *cw_im, double *jw, 
			fftw_complex *gt,double *w_m0_sh, float ***sideband, float ***zeroline, 
			double *spek_ab, double *sumspek_ab, 
			double *spek_sb, double **arrayspek_sb, double *sumspek_sb, 
			double *spek_zl, double **arrayspek_zl, double *sumspek_zl, 
			double *spek_cd, double *sumspek_cd, double *spek_ld, double *sumspek_ld, 
			double *spek_cdmk, double *sumspek_cdmk, double *spek_ldmk, double *sumspek_ldmk, 
			double *sym_ax, fftw_complex *in, fftw_complex *out, fftw_plan plan)
{

     
  matrix2dia(nz, site, sigma, seed, vab, hamilt);       // Hamilton-Matrix aufstellen 

  rs(nz, hamilt, eigval, 1, ev_h);                      // Eispack-Rotine: EW und EV berechnen  
  ewev(nz, eigval, ev_h, eigvec);                       // konvert EV(1D) => EV(2D)
  
  linAbsKoeff(nz, qy, eigvec, alpha, dipst);                //Berechnung von alpha=|dipst|^2, [alpha]=DEBYE^2
  rotatStrength(nz, qy, eigvec, rotst, dipst, rij, delta);
  linDikroKoeff(nz, sym_ax, alpha, qy, ldcff, eigvec, dipst);


  gammaTau(nz, rij, eigval, eigvec, gam, tau, cw_re, jw);   // cw_re wird nicht benutzt!

  omegaShift(nz, eigval, gam, w_m0_sh, cw_im, jw);  

  dynTheoryArray(nz, anz, sideband, zeroline, gt, gam, tau, w_m0_sh, spek_ab, alpha, 
		 spek_sb, arrayspek_sb, spek_zl, arrayspek_zl, 
		 spek_cd, rotst, spek_ld, ldcff, spek_cdmk,  spek_ldmk); // OD + CD + LD spektren
  
  spekSum(anz, spek_ab, sumspek_ab);
  spekSum(anz, spek_sb, sumspek_sb);
  spekSum(anz, spek_zl, sumspek_zl);
  spekSum(anz, spek_cd, sumspek_cd);
  spekSum(anz, spek_ld, sumspek_ld);
  spekSum(anz, spek_cdmk, sumspek_cdmk);
  spekSum(anz, spek_ldmk, sumspek_ldmk);


}

  // pigmExcDistrib(nz, w_m0_sh, npig, nexc, eigvec);


//============================================================================================================================


void calcLinSpecMarkov(int nz, int i, int anz, double *site, double *sigma, int *seed, double **vab, double *hamilt, 
			double *eigval, double **eigvec, double *ev_h, double **qy, double *alpha, double *rotst, double *ldcff, 
			double **delta, double *dipst, double **rij, double **gam, double *tau, 
			double *cw_re, double *cw_im, double *jw, fftw_complex *gt, double *w_m0_sh, 
			double *spek_ab, double *sumspek_ab, double *spek_cd, double *sumspek_cd, double *spek_ld, double *sumspek_ld, 
			double *sym_ax)
{

     
  matrix2dia(nz, site, sigma, seed, vab, hamilt);       // Hamilton-Matrix aufstellen 

  rs(nz, hamilt, eigval, 1, ev_h);                      // Eispack-Rotine: EW und EV berechnen  
  ewev(nz, eigval, ev_h, eigvec);                       // konvert EV(1D) => EV(2D)
  
  linAbsKoeff(nz, qy, eigvec, alpha, dipst);                //Berechnung von alpha=|dipst|^2, [alpha]=DEBYE^2
  rotatStrength(nz, qy, eigvec, rotst, dipst, rij, delta);
  linDikroKoeff(nz, sym_ax, alpha, qy, ldcff, eigvec, dipst);


  gammaTau(nz, rij, eigval, eigvec, gam, tau, cw_re, jw); 

  omegaShift(nz, eigval, gam, w_m0_sh, cw_im, jw);   

  //omegaShiftDia(nz, eigval, gam, w_m0_sh);
  //omegaUnShifted(nz, eigval, w_m0_sh);
  
 
  markovApprox(nz, anz, alpha, rotst, ldcff, tau, w_m0_sh, spek_ab, spek_cd, spek_ld); 


  spekSum(anz, spek_ab, sumspek_ab);
  spekSum(anz, spek_cd, sumspek_cd);
  spekSum(anz, spek_ld, sumspek_ld);


}


//==============================================================================================================


void calcHoleBurnSpectrum(int nz, int pigm_no, int anz, double wstep, double *site, double *sigma, int *seed, 
			  double **vab, double **vab_hb, double *hamilt, double *eigval, 
			  double **eigvec, double *eigval_hb, 
			  double **eigvec_hb, double *ev_h, double **qy, double *alpha, double *dipst, 
			  double **rij, double **gam, double *tau, double *cw_re, double *cw_im, double *jw, 
			  fftw_complex *gt, double *w_m0_sh, double *w_hb_sh, float ***sideband, float ***zeroline, 
			  double **arrayspek_sb, double **arrayspek_zl, 
			  double *spek_ab, double *spek_sb, double *spek_zl, double *spek_hb, double *spek_hb_res, 
			  double *sumspek_hb, double *sumspek_hbab, fftw_complex *in, fftw_complex *out, fftw_plan plan)
{   

  // Hole-Burning-Spektrum: Brennen von Pigment pigm_no und berechnen des Anteils am HB Spektrum

  int k=0, j=0;
  double p_pigm, p_burn, p_burn_res;


  // Calculate resonant HB-Spektrum:

  matrix2diaHBres(nz, site, sigma, seed, w_m0_sh, vab, vab_hb, hamilt, pigm_no); // site energies resonant burning 
  rs(nz, hamilt, eigval_hb, 1, ev_h);                                                // Eispack-Rotine: EW und EV berechnen  
  ewev(nz, eigval_hb, ev_h, eigvec_hb);                                              // konvert EV(1D) => EV(2D)
  
  linAbsKoeff(nz, qy, eigvec_hb, alpha, dipst);                                      //Berechnung des linearen Absorptionskoeff. 
  gammaTau(nz, rij, eigval_hb, eigvec_hb, gam, tau, cw_re, jw); 
  omegaShift(nz, eigval_hb, gam, w_hb_sh, cw_im, jw); 
 
  dynTheoryArrayHB(nz, anz, sideband, zeroline, gt, gam, tau, w_hb_sh, spek_hb_res, alpha, vab, vab_hb);  //resonant
  //  spekDiff(anz, spek_ab, spek_hb_res);   // homogen HB-DiffSpek

  //-------------------------------------------------------------------------------------------------------------------

  // Calculate non-resonant HB-Spektrum:

  matrix2diaHB(nz, site, sigma, seed, w_m0_sh, vab, vab_hb, hamilt, pigm_no); // site energies non-resonant burning 
  rs(nz, hamilt, eigval_hb, 1, ev_h);                                         // Eispack-Rotine: EW und EV berechnen  
  ewev(nz, eigval_hb, ev_h, eigvec_hb);                                       // konvert EV(1D) => EV(2D)

  linAbsKoeff(nz, qy, eigvec_hb, alpha, dipst);                               // Berechnung des linearen Absorptionskoeff. 
  gammaTau(nz, rij, eigval_hb, eigvec_hb, gam, tau, cw_re, jw); 
  omegaShift(nz, eigval_hb, gam, w_hb_sh, cw_im, jw); 

  dynTheoryArrayHB(nz, anz, sideband, zeroline, gt, gam, tau, w_hb_sh, spek_hb, alpha, vab, vab_hb);     //non-resonant
  //  spekDiff(anz, spek_ab, spek_hb);        // homogen HB-DiffSpek

  //--------------------------------------------------------------------------------------------------------------------

  p_pigm = SQ(eigvec[pigm_no][0]);   // Pigm, Exciton 1, Wahrscheinlichkeit P(i,k) = |C(i,k,M=1)|^2
                                     // dass Pigment pigm_no im niedrigsten Exzitonenzustand ist

  // Brennwahrscheinlichkeit P_burn bei der Brennfrequenz
  // E_burn = SPEKMIN + k*wstep  =>  k = (E_burn - SPEKMIN)/wstep;

  k = (int) (E_BURN - SPEKMIN)/wstep;
  
  p_burn_res = arrayspek_zl[0][k];                                           // resonant burning in 0-0-line 

  p_burn = 0.;
  for(j=1; j<nz; j++)
    p_burn += (arrayspek_zl[j][k] + arrayspek_sb[j][k]);
  p_burn += arrayspek_sb[0][k];



  for(j = 0; j<anz; j++)  //Berechnung des HB Spektrums fuer Pigment pig_no
  {
    sumspek_hb[j]   += p_pigm*( p_burn_res*spek_hb_res[j] + p_burn*spek_hb[j] )/(1.0*N); //post-burn
    sumspek_hbab[j] += p_pigm*( p_burn_res + p_burn )*spek_ab[j]/(1.0*N);                //pre-burn
  }
 

}



  // p_burn    = arrayspek_sb[0][k] + arrayspek_zl[1][k] + arrayspek_sb[1][k]; // non-resonant burning 
  // p_burn: Achtung, für mehr als 2 Pigm anpassen! 

  // arrayspek_zl[M][k]                      = |mu_01|^2*ZL_1(w_exc)
  // arrayspek_zl[M][k] + arrayspek_sb[M][k] = |mu_0M|^2*D_M(w_exc)


  //  sumKoeff(i, nz, eigvec, sumkoeff);


