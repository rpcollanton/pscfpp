BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pspg/include.mk
include $(SRC_DIR)/pspg/tests/sweep/sources.mk

TEST=pspg/tests/sweep/Test

all: $(pspg_tests_sweep_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

libs:
	@echo $(LIBS)

run: $(pspg_tests_sweep_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspg/tests/ 

run-quiet: $(pspg_tests_sweep_OBJS) $(BLD_DIR)/$(TEST)
	@$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspg/tests/  > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pspg/tests/sweep/log" > count
	@cat count

clean-outputs:
	rm -f log count 
	find out -type f \( ! -iname ".*" \) -exec rm -f {} \;

clean:
	rm -f $(pspg_tests_sweep_OBJS) $(pspg_tests_sweep_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	$(MAKE) clean-outputs

-include $(pspg_tests_sweep_OBJS:.o=.d)
-include $(pspg_tests_sweep_OBJS:.o=.d)
