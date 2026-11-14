# nug30 with the node trace: recipes/nug30.sh (the acceptance run) plus --trace, so that every
# node's line (QAP_NODE_TRACE, latest/bnb.C) goes to the store with its unit and the frontier's to
# runs/RUN/frontier/: the first tree the viewer can show whole (about 3,700 nodes). Sourced from
# the repository root, as node.sh and the local runs do.
. infra/phase1/recipes/nug30.sh
COORD="$COORD --trace"
