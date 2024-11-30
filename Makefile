# FC= /usr/local/intel/fc/9.0/bin/ifort
FC = ifort
FLAGS= -c  -w   -r8  -O3 

all:    
	(make -f Compile_SAC    FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) ;\
	(make -f Compile_SAC_p  FC="$(FC)" FLAGS="$(FLAGS)"  LIBS="$(LIBS)" ) 

clean: 	
	(make -f Compile_SAC  clean );\
	(make -f Compile_SAC_p  clean );\
	rm *.mod
