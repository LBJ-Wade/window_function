#
# Two-point correlation function calculator
#

CXX       = g++ -std=c++11
EXEC      = rr_brute_force2
DIRS      = $(GSL_DIR)
VPATH     = ../common

all: $(EXEC)

INCLDIRS  = $(foreach dir, $(DIRS), -I$(dir)/include)
LIBDIRS   = $(foreach dir, $(DIRS), -L$(dir)/lib) 
CXXFLAGS := $(INCLDIRS) -I../common

OBJS1 = rr_brute_force2.o read.o histogram2.o
LIBS1 := $(LIBDIRS) -lm -lboost_program_options
LIBS1 += -lgsl -lgslcblas

rr_brute_force2: $(OBJS1)
	$(CXX) $(OBJS1) $(LIBS1) -o $@

rr_brute_force2.o: rr_brute_force2.cpp read.h particle.h histogram2.h
read.o: read.cpp read.h particle.h


.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS1)
