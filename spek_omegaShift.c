/**********************************************************************************************

  Instead of calculating omegaShift discrete for each monte-carlo-loop, we calculate 
  the C_im(w) before starting the MC-loop. 
  The omega-corresponding array number is calculated and the value for C_im(w) out is taken 
  of the array calculated in calcCimOmega. This saves a lot of time. 

  w_m0_sh = w_m0 - gam_MM*e_lamb/hbar + sum_(M!=M) gam_MK*C_im(w_MK)    // Reng&Marc 2002, 4.19 
  

 **********************************************************************************************/

#include "spek_head.h" 

void omegaShift(int nz, double *eigval, double **gamma, double *w_m0_sh, double *cw_im, double *jw){
 
  int j=0, k=0, m=0, n=0;
  
  double wstep = WSTEP;           // Schrittweite [fs-1]
  double sum[nz], wmk[nz][nz];

  for(k=0; k<nz; k++) 
    for(m=0; m<nz; m++) 
      wmk[k][m]=0.0; 
      
  for(j=0; j<nz; j++) 
    sum[j] = 0.0;

  for(m=0; m<nz; m++)
    for(k=0; k<nz; k++)
      wmk[m][k] = (eigval[m]-eigval[k])*WZ2WFS;   // [fs-1] 
    

  for(j=0; j<nz; j++){
    for(k=0; k<nz; k++){
      if (k != j)
      {                                    
	n = NT/2 + round( wmk[j][k]/wstep );    // = w_MK/Schrittweite, + NT/2 wg Verschiebung symmetrisch zu Null 

	sum[j] += gamma[j][k]*cw_im[n];     // [sum] = cm-1
      }
    }
  }

 
  //printf("\n"); 

  // Reorganisation Energie  E_LAMB [cm-1],  w_M0 = eigval   [cm-1]
  
  for(j=0; j<nz; j++)         
  {  
    w_m0_sh[j] = eigval[j] - gamma[j][j]*E_LAMB + sum[j];    // [cm-1] 

    //printf("OmegaShift: %5.0lf cm-1,  Eigenvalue: %8.3lf cm-1, gam*E_lam: %5.1lf, Sum: %5.4lf\n", 
    //w_m0_sh[j], eigval[j], -gamma[j][j]*E_LAMB, sum[j]);
  }
  //printf("\n");

}


//===========================================================================


void omegaShiftDia(int nz, double *eigval, double **gamma, double *w_m0_sh){
 
  int j=0;

  // Shifted Peak Position in DiagonalPart-Approximation
  // Reorganisation Energie  E_lamb /approx  102 1/cm

  for(j=0; j<nz; j++){       
    w_m0_sh[j] = eigval[j] - gamma[j][j]*E_LAMB; 
   }
     
}
//===========================================================================


void omegaUnShifted(int nz, double *eigval, double *w_m0_sh){
 
  int j=0;

  for(j=0; j<nz; j++){       
    w_m0_sh[j] = eigval[j];  
   }
     
}

//==============================================================================

  
