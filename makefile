#Compilation instructions:
# By default, `make` will use `gcc` and the bundled libraries
# If you want to compile against a different copy of the INTEL libraries, adjust MKL_LIB_DIR and INTEL_OMP_DIR
# For example:
# 	make cc=ICC MKL_LIB_DIR=/opt/intel/compilers_and_libraries/linux/mkl/intel64 INTEL_OMP_DIR=/opt/intel/compilers_and_libraries/linux/compiler/intel64

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
MKL_LIB_DIR:=/anaconda/lib
INTEL_OMP_DIR:=/anaconda/lib

# Default libraries and flags
FLAGS:=-DMKL_ILP64 -std=c99
RPATH:=-rpath,${MKL_LIB_DIR},-rpath,${INTEL_OMP_DIR}
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
ifeq (${DEBUG},1)
FLAGS+=-DDEBUG -g -Wall
endif

.PHONY: all clean cython_extension clean_c clean_cython

all: wfes cython_extension

clean: clean_c clean_cython

wfes: ${MODULES} ${HEADERS}
	${CC} ${MODULES} -o $@ -I${INC_DIR} -L${MKL_LIB_DIR} -Wl,${LINKER_FLAGS} ${FLAGS} ${LIBS}

cython_extension: wfes.pyx setup.py
	python setup.py build_ext --inplace

params.txt: generate_params.py
	python $< > $@

clean_c:
	${RM} wfes *.so __pycache__

clean_cython:
	${RM} wfes.c build/ *.so
