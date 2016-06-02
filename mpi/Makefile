CC ?= gcc
MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile) -std=c99
MPI_LINK_FLAGS = $(shell mpicc --showme:link)
INCLUDE = -I../include
LIBS = -fopenmp
DOCS=../docs

PROGS = sweep-tt-multistart

all: $(PROGS)

$(PROGS) : % : %.c
	$(CC) $(MPI_COMPILE_FLAGS) $(LIBS) $(INCLUDE) $< $(MPI_LINK_FLAGS) -o $@

job: sweep-tt-multistart
	qsub -q generic run.sh

clean:
	$(RM) run.sh.* output $(PROGS) output.tt
