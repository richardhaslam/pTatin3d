
#
#  Modified makefile for spmA 
#  Here we enforce that the petsc compilers be used
#


# Include the petsc defined variables, this includes compilers, linkers and compilation flags
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


# Don't use plain gcc, use the compiler defined via petsc
#SPMA_CC = gcc
SPMA_CC = ${PCC}

SPMA_CFLAGS = -g -O0 -Wall

all:
	${SPMA_CC} -c ${SPMA_CFLAGS} spmA.c

driver: spmA.o
	${SPMA_CC} -c ${SPMA_CFLAGS} spmA_driver.c
	${SPMA_CC} ${SPMA_CFLAGS} -o spmA_driver spmA_driver.o spmA.o

clean_all:
	rm -rf spmA.o
	rm -rf spmA_driver.o
	rm -rf spmA_driver

