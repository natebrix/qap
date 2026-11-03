# make -C latest -f makefile -f ../tests/dnn_pool.mk poolcheck
# dnn_pool.C includes dnncut.C for exact objective checks; do not link dnncut.o twice.
POOL_TEST ?= /tmp/qap-dnn-pool-test
.PHONY: poolcheck
poolcheck: $(OBJS) dnn.o
	$(CPP) $(CFLAGS) -I. $(INCLUDES) ../tests/dnn_pool.C $(OBJS) dnn.o -o $(POOL_TEST) $(LDLIBS) $(DNN_LAPACK)
	$(POOL_TEST)
