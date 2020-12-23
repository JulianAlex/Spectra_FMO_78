# Spectra_FMO_78
Calculate linear Spectra of FMO (7 and 8 pigments)

Julian Adolphs 

Programm zu Berechnung der inhomogenen linearen Spektren (Absorption OD, circular Dichroism CD, linear Dichroism LD) 
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
