
void calcCt(fftw_complex *ct, double *jw, fftw_complex *in, fftw_complex *out, fftw_plan plan);
void calcCw(fftw_complex *ct, fftw_complex *cw, double *cw_re, double *cw_im, fftw_complex *in, fftw_complex *out, fftw_plan plan);
void calcGt(fftw_complex *gt, double *jw, fftw_complex *in, fftw_complex *out, fftw_plan plan);
void calcJw(double *jw);
void readCoordinates(int *num, double *x, double *y, double *z, int n, char *in_dat); 
void readTRImerCoord(int *trim_num, double **trim_r, int trim_nz, char *in_dat);
void readCouplings(int nv, int nz, double **vab, char *in_dat);
int countInput(char *);   

void maxnorm(int n, double *ab);
void normSB(double *ab,double *hb);
void koord(int nz, double *x, double *y, double *z, double **na, double **nb, double **nc, double **nd);
void dipolmom(int nz, double **na, double **nb, double **nc, double **nd, double **qx, double **qy);
void pigmCenters(int nz, double **mg, double **nb, double **nd, double **delta, double **rij);
void init(int nz, double *site, double *sigma, double *dipst);
void parameterFit(double *site, int nz, int i, int k, double **parent);
void parameter(double *site, double *sigma, double *dipst, int bclnum, int nz);
void interaction(int nz, double **qy, double **delta, double **vab,double **rij, double *dipst);
void calcDipoleStrength(int n1, double *x1, double *y1, double *z1,double *charge);
void calcDipoleMoment(int nz, int n1, double *x7, double *y7, double *z7,double *charge, double **qy);
void coulombCouplings(int nz, int n1, double *x7, double *y7, double *z7,double *charge, double **vab);
void matrix2dia(int nz, double *site, double *sigma, int *seed, double **vab,double *hamilt);
void matrix2diaHB(int nz, double *site, double *sigma, int *seed, double *w_m0_sh, double **vab, double **vab_hb, 
		  double *hamilt_hb, int pigm_no);
void matrix2diaHBres(int nz, double *site, double *sigma, int *seed, double *w_m0_sh, double **vab, double **vab_hb, 
		     double *hamilt_hb, int pigm_no);
void ewev(int nz, double *eigval, double *ev_h, double **eigvec);
void linAbsKoeff(int nz, double **qy, double **eigvec, double *alpha, double *dipst);
void linDikroKoeff(int nz, double *sym_ax, double *alpha, double **qy, double *ld_koeff, double **eigvec, double *dipst);
void rotatStrength(int nz, double **qy, double **eigvec, double *rotst, double *dipst, double **rij, double **delta);
void trimSymAx(int nz, double **trim_r, double *sym_ax);
void sumKoeff(int i, int nz, double **eigvec, double **sumkoeff);
void spekDiffHB(int anz, double *spek_ab, double *spek_hb, double *spek_hbdif);
void spekSum(int anz, double *spek, double *sumspek);
void reset(double *array, int n);
void pigmExcDistrib(int nz, double *w_m0_sh, double **npig, double **nexc, double **nr_ev);

double mutElement(int mut_elem, int nz, int i, double **offspring, double mutmax);
int selectElement(int nz);
void selectNinsrt(int *irank, int i, int k, int *num_sel );
void notinsert(double *fitness,  int *irank, int i, int k, int nz, int *num_sel, double **parent, double **offspring);
void twoSelect(double *vertlgsfkt, int *irank, int i, int k, int *indx_sel_1,int *num_sel_1, int *indx_sel_2, int *num_sel_2);
void rekomb2(double *fitness, double *vertlgsfkt, int *irank, int i, int k,int nz, int *indx_sel_1, int *num_sel_1, int *indx_sel_2, 
	     int *num_sel_2, double **parent, double **offspring);
void selection(double *vertlgsfkt, int *irank, int i, int k, int *indx_sel, int *num_sel );
void mutation(double *fitness, double *vertlgsfkt, int *irank, int i, int k, int nz, int *indx_sel, 
	      int *num_sel, double **parent, double **offspring);
void rekombination(double *fitness, double *vertlgsfkt, int *irank, int i, int k, int nz, int *indx_sel_1, int *num_sel_1, 
		   int *indx_sel_2, int *num_sel_2, double **parent, double **offspring);
void reproduktion(double *fitness, double *vertlgsfkt, int *irank, int i, int k, int nz, int *indx_sel, int *num_sel, double **parent, double **offspring);
void repromut(double *fitness, double *vertlgsfkt, int *irank, int i, int k, int nz, int *indx_sel, int *num_sel, double **parent, double **offspring);
void geneticOperate(int n1, int n2, int n3, int n4, int nz, int i, int k, int *irank, int *num_sel, int *num_sel_1, int *num_sel_2, 
		    int *indx_sel, int *indx_sel_1, int *indx_sel_2, double *fitness, double *vertlgsfkt, double **parent, double **offspring);

void diffExpSim(int i, int anz, double wstep, double *fitness, double *nr_fitness,double *sumspek_ab, double *exp_ab);
void verteilungsfkt(double *fitness, double *nr_fitness, double *rank_fitn, double *vertlgsfkt, 
		    double *prob, int *irank, int *nr_irank, double **parent);

void gammaTau(int nz, double **rij, double *eigval, double **eigvec, double **gamma, double *tau, double *cw_re, double *jw);
void omegaShift(int nz, double *eigval, double **gamma, double *w_m0_sh, double *cw_im, double *jw);


void dynTheoryArray(int nz, int anz, float ***sideband, float ***zeroline, fftw_complex *gt, double **gam, double *tau, double *w_m0_sh,
		    double *spek_ab, double *alpha, double *spek_sb, double **arrayspek_sb, double *spek_zl, double **arrayspek_zl,  
		    double *spek_cd, double *rotst, double *spek_ld, double *ldcff, double *spek_cdmk, double *spek_ldmk);
void dynTheoryArrayHB(int nz, int anz, float ***sideband, float ***zeroline,
		      fftw_complex *gt, double **gam, double *tau, double *w_m0_sh,
		      double *spek_hb, double *alpha, double **vab, double **vab_hb);

void indexxing(int n, double *arr, int *indx);
void ranking(int n, int *indx, int *irank);
void sorting(int n, double arr[]);

void lineshapeArray(fftw_complex *gt, float ***sideband, float ***zeroline, fftw_complex *in, fftw_complex *out, fftw_plan plan);

void eingabeExpSpek(char *in_dat, double *exp_ab, double *exp_cd, double *exp_ld, int nx);

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
			double *sym_ax, fftw_complex *in, fftw_complex *out, fftw_plan plan); 

void markovApprox(int nz, int anz, double *alpha, double *rotst, double *ld_koeff, 
		  double *tau, double *w_m0_sh, double *nmark_od, double *nmark_cd, double *nmark_ld);

int outputSpectrum(double *sumspek, int anz, char *out_dat);
 
