
pspg_sweep_= \
  pspg/sweep/BasisFieldState.cu \
  pspg/sweep/Sweep.cu \
  pspg/sweep/LinearSweep.cu \
  pspg/sweep/SweepFactory.cu

pspg_sweep_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_sweep_))
pspg_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_sweep_:.cu=.o))

