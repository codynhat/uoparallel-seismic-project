CC ?= gcc
CUDA_CC ?= nvcc
MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile) -std=c99
MPI_LINK_FLAGS = $(shell mpicc --showme:link)
CUDA_LINK_FLAGS = -L/usr/local/packages/cuda/7.5/lib64 -lcudart -lcuda -lstdc++
INCLUDE = -I../include
LIBS = -fopenmp
DOCS=../docs

PROGS = sweep-mpi-omp
CUDA_PROGS = sweep-mpi-cuda

all: $(PROGS)

$(PROGS) : % : %.c
	$(CC) $(MPI_COMPILE_FLAGS) $(LIBS) $(INCLUDE) $< $(MPI_LINK_FLAGS) -o $@

mpi-cuda : % : sweep-mpi-cuda.c sweep-cuda.o
	$(CC) $(MPI_COMPILE_FLAGS) $(LIBS) $(INCLUDE) $< sweep-cuda.o $(CUDA_LINK_FLAGS) $(MPI_LINK_FLAGS) -o sweep-mpi-cuda

sweep-cuda.o : sweep-cuda.cu
	$(CUDA_CC) $(INCLUDE) -c $< -o $@

job: sweep-mpi-omp
	qsub -q generic run.sh

clean:
	$(RM) run.sh.* output $(PROGS) output.tt $(CUDA_PROGS)
