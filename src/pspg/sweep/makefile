#--------------------------------------------------------------------
# Include makefiles

SRC_DIR_REL =../..
include $(SRC_DIR_REL)/config.mk
include $(SRC_DIR)/pspg/include.mk

#--------------------------------------------------------------------
# Main targets 

all: $(pspg_sweep_OBJS) 

includes:
	echo $(INCLUDES)

clean:
	rm -f $(pspg_sweep_OBJS) $(pspg_sweep_OBJS:.o=.d) 

#--------------------------------------------------------------------
# Include dependency files

-include $(pspg_OBJS:.o=.d)
