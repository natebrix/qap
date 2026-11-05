# Build the solver and run the checks. See RUNNING.md and CONTRIBUTING.md.
#
#   make               build src/MESCHACH/meschach.a and latest/bnbdrv
#   make test          brute-force, resume and parallel tests, then the regression check
#   make bruteforce    solver vs full enumeration on small random instances (~10 s)
#   make regression    quick suite vs bench/reference/ (~15 s)
#   make resume        stop/resume a run in slices; totals must not change (~10 s)
#   make parallel      bench/parallel.py matches the serial solver (~20 s)
#   make dnntest       latest/dnntest, the DNN bound's driver (needs LAPACK, see latest/makefile)
#   make dnncheck      validate the DNN bound (tests/dnn_test.py, ~2 min)
#   make continuationcheck  same-node DNN continuation and exact-bound checks
#   make clean

PYTHON ?= python3

.PHONY: all solver test bruteforce resume parallel coord regression clean dnntest dnncheck shipcheck

all: solver

solver: src/MESCHACH/meschach.a
	$(MAKE) -C latest

# rebuilt from scratch when its makefile changes (its .o files depend only on the .c files)
src/MESCHACH/meschach.a: src/MESCHACH/makefile
	$(MAKE) -C src/MESCHACH cleanup
	$(MAKE) -C src/MESCHACH basic

test: bruteforce resume parallel coord regression

bruteforce: solver
	$(PYTHON) tests/brute_force.py

resume: solver
	$(PYTHON) tests/resume_test.py

parallel: solver
	$(PYTHON) tests/parallel_test.py

# the coordinator (coord/): totals equal serial under worker kills, a
# coordinator restart and a drain; the AWS scenario runs when boto3 and moto
# are importable (else it is skipped)
coord: solver
	$(PYTHON) tests/coord_test.py

regression: solver
	$(PYTHON) bench/bench.py check

dnntest: src/MESCHACH/meschach.a
	$(MAKE) -C latest dnntest

dnncheck: dnntest
	$(PYTHON) tests/dnn_test.py

# the parent's DNN state through a stop and resume (state v6, QAP_DNN_SHIP), about 4 minutes
shipcheck: solver
	$(PYTHON) tests/ship_test.py

clean:
	$(MAKE) -C latest clean
	rm -f latest/bnbdrv

.PHONY: continuationcheck
continuationcheck: src/MESCHACH/meschach.a
	$(MAKE) -C latest dnncontinuation
	./latest/dnncontinuation
