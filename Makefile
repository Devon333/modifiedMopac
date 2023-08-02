##############################################################################
#
#   Compile and build MOPAC2016 for Linux   
#
##############################################################################
#
#  The following set of definitions apply to this OS.
#
##############################################################################
MKLROOT= /share/software/languages/INTEL-2015/composer_xe_2015
ifort = $(MKLROOT)/bin/ifort
MATH_LINK       = -Wl,--start-group  $(MKLROOT)/mkl/lib/intel64/libmkl_intel_lp64.a \
                  $(MKLROOT)/mkl/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/mkl/lib/intel64/libmkl_core.a \
                  -Wl,--end-group -openmp -lpthread -lm
FORTRAN_LINKER  = $(MKLROOT)/bin/ifort
FORTRAN_COMPILER= $(MKLROOT)/bin/ifort -o $@ -c -lpthread -lstdc++ -O3 -I/$(OBJ) -module $(OBJ) -diag-disable 8290 -diag-disable 8291

MOPAC_SRC        = ./src
OBJ              = ./obj
O                = o 
EXE              = ./MOPAC2016.exe
include MOPAC_Makefile_files.txt
