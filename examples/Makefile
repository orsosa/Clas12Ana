.DELETE_ON_ERROR:

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

HIPOCFLAGS  := -I$(HIPO4INC)
HIPOLIBS    := -L$(HIPO4LIB) -lhipo4

LZ4LIBS     := -L$(LZ4DIR)/lib -llz4
LZ4INCLUDES := -I$(LZ4DIR)/lib

SLIB_DIR  := -L$(CLAS12ANA_H4)/shlib -lTIdentificatorCLAS12.so
STLIB_DIR := -L$(CLAS12ANA_H4)/stlib -lTIdentificatorCLAS12.a
INC_DIR   := -I$(CLAS12ANA_H4)/include

CXX       := g++
CXXFLAGS  += -Wall -fPIC $(ROOTCFLAGS) -std=c++11
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)
SOFLAGS   = -Wl,-soname,$(notdir $@) -shared

AR	  = ar
ARFLAGS	  = -cvr #create,verbose,quick (don't check for replacement, otherwise use r instead)


all:  getSchema get_simple_tuple particle_mix

%: %.o
	@echo  $(ANALYSERLIBS)
	$(CXX) -o $@ $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS) $(ANALYSERLIBS)

%.o: %.cxx
	$(CXX) -c -o $@ $< $(ROOTCFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES) $(INC_DIR)

