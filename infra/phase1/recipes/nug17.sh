# nug17, the smoke test: QPB only (no DNN), about 20 CPU s. One-second units
# push a few dozen of them through S3, SQS and DynamoDB, so a first launch
# checks the whole path for cents before nug30.
INSTANCE=qaplib/nug17.dat
PARAMS=PARAM/default.par
OPT=1732          # the known optimum, which the check expects
UB=1732           # the incumbent given to the solver (--no-ub drops it)
EXTRA=""
ENV=""
COORD="--budget 1 --target-unit-cpu 1 --steal-after 2 --checkpoint-every 0.5"
