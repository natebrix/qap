CPP = g++
BNB = bnbdrv

#CF = -pg -O4
#CF = -g -O4
# -march=native lets the compiler use this CPU's instruction set (~2%).
# -ffp-contract=off keeps it from fusing multiply-adds, so the arithmetic, and
# therefore node counts, stay identical across machines (see RUNNING.md).
CF = -O3 -march=native -ffp-contract=off

CFLAGS = $(CF) -DBNB_STANDALONE

# The DNN bound (dnn.C; bnbdrv -x) needs LAPACK, found as described at
# the dnntest target below. Build without it: make DNN=0.
DNN ?= 1
ifeq ($(DNN),1)
CFLAGS += -DHAVE_DNN
DNN_OBJ = dnn.o
DNN_LIBS = $(DNN_LAPACK)
endif

MESCHACH = ../src/MESCHACH
MESCHACHLIB = ../src/MESCHACH/meschach.a

INCLUDES = -I$(MESCHACH) -I$(EVB3) -I./
 
LDLIBS = -lm $(MESCHACHLIB)

OBJS = qpbnd.o qpbnd2.o lap.o dlap.o qap.o assign.o link1.o \
       solution.o branch.o bnb.o util.o linalg.o fw.o estimate.o svd.o \
       sa_qap.o stats.o console.o lgb_train.o state.o

all:    bnb

clean:
	rm -f *.o *~ a.out dnntest

.C.o:	
	$(CPP) $(CFLAGS) $(INCLUDES) -c $?

bnb:	$(OBJS) $(DNN_OBJ) driver.o
	$(CPP) $(CFLAGS) $(OBJS) $(DNN_OBJ) driver.o -o $(BNB) $(LDLIBS) $(DNN_LIBS)


# ---- DNN bound (dnn.C, dnn.h) and its test driver: make dnntest ----
# dnn.C needs LAPACK/BLAS. If the scipy-openblas32 wheel
# is installed (pip install scipy-openblas32), its OpenBLAS is used: 32-bit
# integers, symbols prefixed scipy_, about 7x faster than the reference
# LAPACK. Otherwise the system liblapack.so.3 / libblas.so.3 are linked.
# Override both with e.g.
#   make dnntest DNN_LAPACK="-lopenblas" DNN_LAPACK_DEFS="-DDNN_OPENBLAS"
# Run make clean after switching: dnn.o depends on DNN_LAPACK_DEFS.
PYTHON ?= python3
OPENBLAS32_DIR = $(shell $(PYTHON) -c "import scipy_openblas32 as o; print(o.get_lib_dir())" 2>/dev/null)
DNN_LAPACK = $(if $(OPENBLAS32_DIR),-L$(OPENBLAS32_DIR) -lscipy_openblas -Wl$(comma)-rpath$(comma)$(OPENBLAS32_DIR),-l:liblapack.so.3 -l:libblas.so.3)
DNN_LAPACK_DEFS = $(if $(OPENBLAS32_DIR),-DDNN_BLAS_PREFIX=scipy_ -DDNN_OPENBLAS,)
comma := ,

dnn.o: dnn.C dnn.h
	$(CPP) $(CFLAGS) $(DNN_LAPACK_DEFS) $(INCLUDES) -c dnn.C

dnntest.o: dnntest.C dnn.h
	$(CPP) $(CFLAGS) $(INCLUDES) -c dnntest.C

dnntest: $(OBJS) dnn.o dnntest.o
	$(CPP) $(CFLAGS) $(OBJS) dnn.o dnntest.o -o dnntest $(LDLIBS) $(DNN_LAPACK)
