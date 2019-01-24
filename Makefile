# Makefile Spektrum
#
#CC_FLAGS = -g -ggdb -Wall -m64 -O3     #-pg -Wall -m64 -O3
#
CC_FLAGS = -g -Wall   

LD_FLAGS =  -lfftw3 -lm

OBJ = spek_main.o spek_correlJw.o spek_dynTheory.o spek_gammaTau.o spek_input.o\
 spek_lineshape.o spek_omegaShift.o spek_output.o spek_spectra.o spek_subworx.o\
 eispack.o normal.o\

SRC = spek_main.c spek_correlJw.c spek_dynTheory.c spek_gammaTau.c spek_input.c\
 spek_lineshape.c spek_omegaShift.c spek_output.c spek_spectra.c spek_subworx.c\
 eispack.c normal.c\

HDR = spek_funk.h spek_head.h
 

CC = gcc     # Compiler

EXECUTABLE = spektrum



$(EXECUTABLE) : $(SRC) $(HDR)                 #$(OBJ)
#
	$(CC) $(CC_FLAGS) -o spektrum $(SRC) $(LD_FLAGS)
#
#
##############################################################################
#
# GPROF:
#
# "gcc -g -Wall -o spektrum -pg -lm spek.c -pg"
#
# compile and run program ( => gprof.out )
#
# run gprof:  "gprof ./spektrum > temp.out"
#

