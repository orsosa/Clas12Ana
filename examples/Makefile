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
INC_DIR   := -I$(CLAS12ANAHIPO4)/include -I.

CXX       := g++
CXXFLAGS  += -Wall -fPIC $(ROOTCFLAGS) -std=c++11 -g -O0
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

AR	  = ar
ARFLAGS	  = -cvr #create,verbose,quick (don't check for replacement, otherwise use r instead)

all: checkdirs slib/libData.so getSchema get_simple_tuple get_simple_tuple_rich get_simple_tuple_th get_simple_tuple_th_mf particle_mix particle_mix_tree getBankInfo dumpTrajBank dumpTrigger get_simple_tree_th_mf

checkdirs: dict slib

dict slib:
	@mkdir -p $@

%: %.o
	$(CXX) -o $@ $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS)  $(ROOTLIBS) $(CLAS12ANALIB) -Lslib -lData

%.o: %.cxx 
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(ROOTCFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES) $(INC_DIR)

dict/datadict.cxx: PARTDATA.h DETDATA.h LinkDef.h 
	rootcling -f $@ -c $(ROOTCFLAGS) $(HIPOCFLAGS) -p $^

slib/libData.so: dict/datadict.cxx PARTDATA.cxx DETDATA.cxx
	g++ -shared -fPIC -o $@ $(ROOTLDFLAGS) $(ROOTCFLAGS) $(HIPOCFLAGS) $(INC_DIR) $(ROOTCFLAGS) $^
	cp dict/datadict_rdict.pcm slib/.

