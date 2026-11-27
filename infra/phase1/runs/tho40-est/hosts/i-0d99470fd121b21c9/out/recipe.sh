# run tho40-est: the recipe launch.sh sent (runs/tho40-est/recipe.sh on the laptop), 128 slots
declare -- INSTANCE="qaplib/tho40t.dat"
declare -- PARAMS="../infra/phase1/recipes/tho40c.par"
declare -- OPT="240516"
declare -- UB="240518"
declare -- EXTRA="-s SYM/tho40t.sym"
declare -- ENV=""
declare -- COORD="--budget 1200 --target-unit-cpu 1200 --steal-after 300 --checkpoint-every 60 --ckpt-upload 300 --trace"
