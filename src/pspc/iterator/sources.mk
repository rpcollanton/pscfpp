pspc_iterator_= \
  pspc/iterator/AmIterator.cpp \
  pspc/iterator/FilmIterator.cpp \
  pspc/iterator/IteratorFactory.cpp 

pspc_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_iterator_:.cpp=.o))

