include ../standard_defs.mk

MEX = $(MATLAB_ROOT_BIN)mex -v CXXFLAGS="\$$CXXFLAGS $(OPENMP)" LDFLAGS="\$$LDFLAGS $(OPENMP)"    
INCL += -I../ -I../rngstream

all: rgig_par

rgig_par_mem_check: rgig_par_mem_check.cpp gig_link.o ../rngstream/RngStream.o ../gig_par.o
	$(CC) -g3 $(CFLAGS) $(INCL) $^ $(LIBS) 

rgig_par: rgig_par.cpp gig_link.o ../rngstream/RngStream.o ../gig_par.o
	$(MEX) $^  $(INCL) $(LIBS)

Egig_par: Egig_par.cpp ../rngstream/RngStream.o ../gig_par.o
	$(MEX) $^ $(INCL)  $(LIBS)
../rngstream/RngStream.o:
	- ( cd ../rngstream/ && $(MAKE))
../gig_par.o:
	- ( cd ../ && $(MAKE))
gig_link.o: gig_link.cpp  
	$(CC)  $(CXXFLAGS) $(INCL) -c $^ -o $@ 
.PHONY: clean
clean:
	rm -f rgig_par *.*~ *.o
	- ( cd ../rngstream/ && $(MAKE) clean)

