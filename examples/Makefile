.DELETE_ON_ERROR:

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

ifndef CLAS12TOOL
$(error CLAS12TOOL is not set)
endif

HIPO4INC := $(if $(HIPO4INC),$(HIPO4INC),$(CLAS12TOOL)/Hipo4 )
HIPO4LIB := $(if $(HIPO4LIB),$(HIPO4LIB),$(CLAS12TOOL)/lib )


HIPOCFLAGS  := -I$(HIPO4INC)
HIPOLIBS    := -L$(HIPO4LIB) -lhipo4

LZ4LIBS     := -L$(LZ4DIR)/lib -llz4
LZ4INCLUDES := -I$(LZ4DIR)/lib

CLAS12ANALIB  := -L$(CLAS12ANAHIPO4)/shlib -lTIdentificatorCLAS12
#CLAS12ANASTLIB := -L$(CLAS12ANAHIPO4)/stlib -lTIdentificatorCLAS12
INC_DIR   := -I$(CLAS12ANAHIPO4)/include

CXX       := g++
CXXFLAGS  += -Wall -fPIC $(ROOTCFLAGS) -std=c++11
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

AR	  = ar
ARFLAGS	  = -cvr #create,verbose,quick (don't check for replacement, otherwise use r instead)


all:  getSchema get_simple_tuple particle_mix

%: %.o
	$(CXX) -o $@ $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS)  $(ROOTLIBS) $(CLAS12ANALIB) 

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(ROOTCFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES) $(INC_DIR)

