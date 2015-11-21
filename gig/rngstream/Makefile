include ../standard_defs.mk
CC = g++
IEXTLIBS =   -fopenmp -I/opt/local/include -I$(CURDIR)/

ifeq ($(shell uname), Linux)
	CXXFLAGS += -std=c++0x
endif

all: RngStream.o 

%.o : %.cpp
	$(CC)  $(CXXFLAGS) $(IEXTLIBS) -c $^ -o $@ 
	
.PHONY: clean
clean:
	rm -f *.o *~ 


