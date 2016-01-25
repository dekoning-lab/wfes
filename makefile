#Compilation instructions:
# By default, `make` will use `gcc` and the bundled libraries
# If you want to compile against a different copy of the INTEL libraries, adjust MKL_LIB_DIR and INTEL_OMP_DIR
# For example:
# 	make cc=ICC MKL_LIB_DIR=/opt/intel/compilers_and_libraries/linux/mkl/intel64 INTEL_OMP_DIR=/opt/intel/compilers_and_libraries/linux/compiler/intel64

CC=gcc
RM=rm -rf
DEBUG=0

SRC_DIR:=src
MODULES:=${SRC_DIR}/*.h
# Paths to MKL libraries (shared by default)
MKL_LIB_DIR:=lib
INTEL_OMP_DIR:=lib

# Default libraries and flags
FLAGS:=-DMKL_ILP64
RPATH:=-rpath,${MKL_LIB_DIR},-rpath,${INTEL_OMP_DIR}
LIBS:=-lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lmkl_avx -lmkl_vml_avx -lm -lrt

# ICC-specific flags
ifeq (${CC},icc)
FLAGS+=-qopenmp -std=c99
LINKER_FLAGS=${RPATH}
endif

# GCC-specific flags
ifeq (${CC},gcc)
FLAGS+=-m64 -fopenmp --std=c99
LINKER_FLAGS=--no-as-needed,${RPATH}
LIBS+=-liomp5 -ldl
endif

# DEBUG build flags
ifeq (${DEBUG},1)
FLAGS+=-DDEBUG -g -Wall
endif

.PHONY: all clean

all: wfes libwfes.so


wfes: ${SRC_DIR}/wfes.c ${MODULES}
	${CC} $< -o $@ -L${MKL_LIB_DIR} -Wl,${LINKER_FLAGS} ${FLAGS} ${LIBS}

libwfes.so: ${SRC_DIR}/wfes.c ${MODULES}
	${CC} -shared -fPIC $< -o $@ -L${MKL_LIB_DIR} -Wl,${LINKER_FLAGS} ${FLAGS} ${LIBS}

clean:
	${RM} wfes libwfes.so __pycache__
