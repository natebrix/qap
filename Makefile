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
#   make clean

PYTHON ?= python3

.PHONY: all solver test bruteforce resume parallel regression clean dnntest dnncheck

all: solver

solver: src/MESCHACH/meschach.a
	$(MAKE) -C latest

src/MESCHACH/meschach.a:
	$(MAKE) -C src/MESCHACH basic

test: bruteforce resume parallel regression

bruteforce: solver
	$(PYTHON) tests/brute_force.py

resume: solver
	$(PYTHON) tests/resume_test.py

parallel: solver
	$(PYTHON) tests/parallel_test.py

regression: solver
	$(PYTHON) bench/bench.py check

dnntest: src/MESCHACH/meschach.a
	$(MAKE) -C latest dnntest

dnncheck: dnntest
	$(PYTHON) tests/dnn_test.py

clean:
	$(MAKE) -C latest clean
	rm -f latest/bnbdrv
