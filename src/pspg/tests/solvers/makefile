BLD_DIR_REL =../../..
include $(BLD_DIR_REL)/config.mk
include $(SRC_DIR)/pspg/include.mk
include $(SRC_DIR)/pspg/tests/solvers/sources.mk

TEST=pspg/tests/solvers/Test

all: $(pspg_tests_solvers_OBJS) $(BLD_DIR)/$(TEST)

includes:
	@echo $(INCLUDES)

libs:
	@echo $(LIBS)

run: $(pspg_tests_solvers_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspg/tests/solvers/

run-quiet: $(pspg_tests_solvers_OBJS) $(BLD_DIR)/$(TEST)
	$(BLD_DIR)/$(TEST) $(SRC_DIR)/pspg/tests/solvers/ > log
	@echo `grep failed log` ", "\
              `grep successful log` "in pspg/tests/solvers/log" > count
	@cat count

clean-outputs:
	rm -f log count 

clean:
	rm -f $(pspg_tests_solvers_OBJS) $(pspg_tests_solvers_OBJS:.o=.d)
	rm -f $(BLD_DIR)/$(TEST) $(BLD_DIR)/$(TEST).d
	$(MAKE) clean-outputs

-include $(pspg_tests_solvers_OBJS:.o=.d)
