pspg_compressor_= \
  pspg/compressor/CompressorFactory.cu \
  pspg/compressor/AmCompressor.cu

pspg_compressor_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_compressor_))
pspg_compressor_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_compressor_:.cu=.o))

