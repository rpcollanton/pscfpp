pspg_tests_field_=pspg/tests/field/Test.ccu

pspg_tests_field_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_tests_field_))
pspg_tests_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_tests_field_:.ccu=.o))
