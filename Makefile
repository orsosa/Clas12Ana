.DELETE_ON_ERROR:

include Makefile_top

SRC_FILES  = $(wildcard $(SRC_DIR)/*.cxx)
SRC_OBJ = $(SRC_FILES:$(SRC_DIR)/%.cxx=$(OBJ_DIR)/%.o)

all: checkdirs $(SLIB_DIR)/$(SH_LIB) app

$(SLIB_DIR)/$(SH_LIB): $(SRC_OBJ)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@	

$(STLIB_DIR)/$(ST_LIB): $(SRC_OBJ)
	$(AR) $(ARFLAGS)  $@  $^
	ranlib $@

app:
	@make -C examples;

%: $(OBJ_DIR)/%.o
	$(CXX) -o $@ $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx
	@echo $(SRC_OBJ)	
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(ROOTCFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES) -I$(INC_DIR)

checkdirs: $(SLIB_DIR) $(OBJ_DIR) $(DICT_DIR) $(DEP_DIR) $(STLIB_DIR)

$(SLIB_DIR) $(OBJ_DIR) $(DICT_DIR) $(DEP_DIR) $(STLIB_DIR):
	@mkdir -p $@

clean:
	@echo 'Removing all build files'
	rm -rf *.o $(OBJ_DIR)/*.o $(SLIB_DIR)/$(SH_LIB)
	@make clean -C examples;

%.o: %.cc
	$(CXX) -c $< -O2 $(ROOTCFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES)
