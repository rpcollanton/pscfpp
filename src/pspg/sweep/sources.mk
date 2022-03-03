
pspg_sweep_= \
  pspg/sweep/FieldState.cpp \
  pspg/sweep/BasisFieldState.cpp \
  pspg/sweep/Sweep.cpp \
  pspg/sweep/LinearSweep.cpp \
  pspg/sweep/SweepFactory.cpp

pspg_sweep_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_sweep_))
pspg_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_sweep_:.cpp=.o))

