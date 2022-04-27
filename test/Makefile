# nourdine.bah@crick.ac.uk

CC=g++
STD=-std=c++14
OPTIM=-O3
FLAGS=-w -pedantic -fopenmp
LIBS=-lz -lrt -lpthread
OPTS=-DNDEBUG
OPTS+=-DSEQAN_HAS_ZLIB=1
OPTS+=-DSEQAN_ENABLE_DEBUG=0
OPTS+=-DSEQAN_ENABLE_TESTING=0
ARGS=$(STD) $(OPTIM) $(FLAGS) $(INC) $(LIBS) $(OPTS)

data.o:
	$(CC) data.cpp $(ARGS) -c -o data.o

create: data.o
	$(CC) data.o create.cpp $(ARGS) -o create

clean:
	rm -fv data.o
	rm -fv create

