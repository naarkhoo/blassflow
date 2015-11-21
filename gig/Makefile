include standard_defs.mk
INCL += -I$(CURDIR)/rngstream



all:  gig_par.o

.PHONY: matlab
matlab:
	- (cd $(CURDIR)/matlab/ && $(MAKE)) 

gig_par.o : gig_par.cpp
	$(CC)  $(CXXFLAGS) $(INCL) $(OPT) -c gig_par.cpp
%.o : %.c
	$(CC)  $(CXXFLAGS) $(INCL) $(OPT) -c $^ -o $@ 


/rngstream/RngStream.o: 
	- ( cd $(CURDIR)/rngstream/ && $(MAKE))


.PHONY: clean
clean:
	rm -f *.o *~
	- ( cd $(CURDIR)/rngstream/ && $(MAKE) clean)
	- (cd $(CURDIR)/matlab/ && $(MAKE) clean)  

