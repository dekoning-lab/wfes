#Compilation instructions:
# If you want to compile against a different copy of the INTEL libraries, adjust MKL_LIB_DIR and INTEL_OMP_DIR
# For example:
# 	make cc=icc MKL_LIB_DIR=/opt/intel/compilers_and_libraries/linux/mkl/intel64 INTEL_OMP_DIR=/opt/intel/compilers_and_libraries/linux/compiler/intel64

uname:=$(shell uname)

ifeq (${uname},Darwin)
CC=clang
else ifeq (${uname},Linux)
CC=gcc
endif

RM=rm -rf
DEBUG=0

SRC_DIR:=src
INC_DIR:=include
MODULES:=${SRC_DIR}/*.c
HEADERS:=${INC_DIR}/*.h

# Paths to MKL libraries (shared by default)
MKL_LIB_DIR:=./lib/${uname}
INTEL_OMP_DIR:=./lib/${uname}

# Path to anaconda MKL libraries (uncomment for building the python extension)
# MKL_LIB_DIR:=${HOME}/anaconda3/lib
# INTEL_OMP_DIR:=${HOME}/anaconda3/lib

# Default libraries and flags
FLAGS:=-DMKL_ILP64 -std=c99
RPATH:=-rpath,${MKL_LIB_DIR},-rpath,'$$ORIGIN/${INTEL_OMP_DIR}'
LIBS:=-lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lmkl_avx -lmkl_vml_avx -lm

# ICC-specific flags
ifeq (${CC},icc)
FLAGS+=-qopenmp -std=c99
LINKER_FLAGS=${RPATH}
endif

# GCC-specific flags
ifeq (${CC},gcc)
FLAGS+=-m64 -fopenmp
LINKER_FLAGS=--no-as-needed,${RPATH}
LIBS+=-liomp5 -ldl -lrt
endif

ifeq (${CC},clang)
FLAGS+=-m64
LIBS+=-liomp5 -ldl
endif

# DEBUG build flags
# use `make DEBUG=1`
ifeq (${DEBUG},1)
FLAGS+=-DDEBUG -g -Wall
endif

# use `make MATRIX_PRINT=1`
ifeq (${MATRIX_PRINT},1)
FLAGS+=-DMATRIX_PRINT
endif

C_RED=\033[0;31m
C_NONE=\033[0m

.PHONY: all clean cython_extension clean_c clean_cython

all: wfes

clean: clean_c clean_cython

wfes: ${MODULES} ${HEADERS}
	${CC} ${MODULES} -o $@ -I${INC_DIR} -L${MKL_LIB_DIR} -Wl,${LINKER_FLAGS} ${FLAGS} ${LIBS}
	@ if [ ${uname} = "Darwin" ]; then echo "${C_RED}[!] please export DYLD_LIBRARY_PATH=${MKL_LIB_DIR}:\$$DYLD_LIBRARY_PATH${C_NONE}"; fi;

cython_extension: wfes.pyx setup.py
	- python setup.py build_ext --inplace

params.txt: generate_params.py
	python $< > $@

clean_c:
	${RM} wfes *.so __pycache__ wfes.dSYM

clean_cython:
	${RM} wfes.c build/ *.so
