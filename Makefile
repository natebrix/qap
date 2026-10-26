# Build the solver and run the checks. See RUNNING.md and CONTRIBUTING.md.
#
#   make               build src/MESCHACH/meschach.a and latest/bnbdrv
#   make test          brute-force tests, resume test, then the regression check
#   make bruteforce    solver vs full enumeration on small random instances (~10 s)
#   make regression    quick suite vs bench/reference/ (~15 s)
#   make resume        stop/resume a run in slices; totals must not change (~10 s)
#   make clean

PYTHON ?= python3

.PHONY: all solver test bruteforce resume regression clean

all: solver

solver: src/MESCHACH/meschach.a
	$(MAKE) -C latest

src/MESCHACH/meschach.a:
	$(MAKE) -C src/MESCHACH basic

test: bruteforce resume regression

bruteforce: solver
	$(PYTHON) tests/brute_force.py

resume: solver
	$(PYTHON) tests/resume_test.py

regression: solver
	$(PYTHON) bench/bench.py check

clean:
	$(MAKE) -C latest clean
	rm -f latest/bnbdrv
