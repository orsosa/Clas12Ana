.DELETE_ON_ERROR:

include Makefile_top

all:  $(SLIB_DIR)/$(SH_LIB) examples

$(SLIB_DIR)/$(SH_LIB): $(SRC_OBJ)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@	

$(STLIB_DIR)/$(ST_LIB): $(SRC_OBJ)
	$(AR) $(ARFLAGS)  $@  $^
	ranlib $@

examples:
	@cd examples;
	@make;

%: $(OBJ_DIR)/%.o
	$(CXX) -o $@ $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx
	$(CXX) -c -o $@ $< $(ROOTCFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES) -I$(INC_DIR)

clean:
	@echo 'Removing all build files'
	@rm -rf *.o getSchema

%.o: %.cc
	$(CXX) -c $< -O2 $(ROOTCFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES)
