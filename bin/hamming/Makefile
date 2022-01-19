
#module load CUDA/10.1.105
#module load GCCcore/8.2.0

CC=g++
INC=-I include
LIBS=-lOpenCL
STD=-std=c++11
#FLAGS=-W -Wall -g
FLAGS=-g -w
OPENCL=-I opencl/include -L opencl/lib

obj/utils.o:
	$(CC) $(INC) $(LIBS) $(STD) $(OPENCL) $(FLAGS) \
		-c \
		-o obj/utils.o \
		lib/utils.c

obj/opencl_kernel.o:
	$(CC) $(INC) $(LIBS) $(STD) $(OPENCL) $(FLAGS) \
		-c \
		-o obj/opencl_kernel.o \
		lib/opencl_kernel.c

obj/opencl_program.o:
	$(CC) $(INC) $(LIBS) $(STD) $(OPENCL) $(FLAGS) \
		-c \
		-o obj/opencl_program.o \
		lib/opencl_program.c

obj/opencl_buffer.o:
	$(CC) $(INC) $(LIBS) $(STD) $(OPENCL) $(FLAGS) \
		-c \
		-o obj/opencl_buffer.o \
		lib/opencl_buffer.c

obj/opencl.o:
	$(CC) $(INC) $(LIBS) $(STD) $(OPENCL) $(FLAGS) \
		-c \
		-o obj/opencl.o \
		lib/opencl.c

obj/methods.o:
	$(CC) $(INC) $(LIBS) $(STD) $(OPENCL) $(FLAGS) \
		-c \
		-o obj/methods.o \
		lib/methods.cpp

hamming: obj/utils.o \
	obj/opencl_kernel.o \
	obj/opencl_program.o \
	obj/opencl_buffer.o \
	obj/opencl.o \
	obj/methods.o
	$(CC) $(INC) $(LIBS) $(STD) $(OPENCL) $(FLAGS) \
		-o hamming \
		obj/utils.o \
		obj/opencl_kernel.o \
		obj/opencl_program.o \
		obj/opencl_buffer.o \
		obj/opencl.o \
		obj/methods.o \
		src/hamming.cpp

