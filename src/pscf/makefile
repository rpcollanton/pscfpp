#------------------------------------------------------------------------
# src/pscf/makefile
#
# The repository copy of this namespace level makefile is in the src/pscf
# directory. This is copied verbatim by the setup script to create the 
# copy in the bld/pscf directory.
#------------------------------------------------------------------------
# Include makefile fragments

SRC_DIR_REL =..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pscf/include.mk

#-----------------------------------------------------------------------
# Main targets 

all: $(pscf_OBJS) $(pscf_LIB)

clean:
	rm -f $(pscf_OBJS) $(pscf_OBJS:.o=.d) $(pscf_LIB)
	cd tests; $(MAKE) clean

veryclean:
	$(MAKE) clean
	-rm -f *.o */*.o
	-rm -f *.d */*.d
	-rm -f lib*.a
ifeq ($(BLD_DIR),$(SRC_DIR))
	-rm -f paths.h
endif


#-----------------------------------------------------------------------
# Include dependency files

-include $(pscf_OBJS:.o=.d)
