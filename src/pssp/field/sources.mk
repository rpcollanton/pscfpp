pssp_field_= \
  pssp/field/RMeshField.cpp \
  pssp/field/KMeshField.cpp \
  pssp/field/FFT.cpp 

pssp_field_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_field_))
pssp_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_field_:.cpp=.o))
