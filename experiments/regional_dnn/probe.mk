# Run from qap/latest with both this file and the repository makefile.
# RESEARCH_DIR is the absolute directory containing region_probe.C.
region_probe: $(OBJS) dnn.o $(RESEARCH_DIR)/region_probe.C
	$(CPP) $(CFLAGS) -I. $(INCLUDES) $(RESEARCH_DIR)/region_probe.C $(OBJS) dnn.o -o $(RESEARCH_DIR)/region_probe $(LDLIBS) $(DNN_LAPACK)
