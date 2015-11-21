CC = /opt/local/bin/g++-mp-4.8  #g++ 
## Use OPENMP?
OPENMP = -fopenmp
MATLAB_ROOT_BIN = 
#/Applications/MATLAB_R2013a.app/bin/



CXXFLAGS = -O3  -frounding-math -fPIC  
ifeq ($(shell uname), Linux)
	INCL     =  -I/usr/include 
	LIBS 	 = -L/usr/lib/
endif

ifeq ($(shell uname), Darwin)
	INCL     =  -I/opt/local/include 
	LIBS 	 = -L/usr/local/lib/ -L/opt/local/lib/
endif
LIBS +=  -lgsl -lm -lblas #gsl is only needed for Egig_par so it can be removed
						  #if Egig is removed	