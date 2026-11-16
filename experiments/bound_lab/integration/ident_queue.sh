#!/bin/bash
# ident_queue.sh BUILD : the switch-off identity runs other than the wil50 slice, one after another for one build
# (the wil50 slice was started by hand: QAP_DNN_POOL_CHECK=20 ./ident.sh wil50_a220 BUILD -p ../experiments/wil50_feasibility/fullrun/fullrun.par
#  -s SYM/wil50t.sym -u 48816 -F 21:22,31:23,40:15 ../qaplib/wil50t.dat, as experiments/estimates/wil50/subtrees/a220)
H=$(cd "$(dirname "$0")" && pwd); B=$1
T="-s SYM/tho40t.sym -u 240516"
# nug20, full search, its own params file, no incumbent given; then with the DNN bound at depths 2-6
$H/ident.sh nug20 $B ../qaplib/nug20.dat
$H/ident.sh nug20_dnn $B -x 2:6:1:300:1 -u 2570 ../qaplib/nug20.dat
# tho40, production settings (no cut phase): the depth-5 subtree "a" of rule2_ab (404 nodes there)
$H/ident.sh tho40_a $B -p ../experiments/estimates/tho40/tho40.par $T -F 26:10,8:0,2:32,18:39,35:24 ../qaplib/tho40t.dat
# tho40 with the grid cut phase in CP form and pools, keys that exist before the merge (tho40_search arm z030; 38 nodes there)
$H/ident.sh tho40_t3_z030 $B -p ../experiments/bound_lab/tho40_search/par/z030.par $T -F 26:10,8:0,38:32 ../qaplib/tho40t.dat
# tai30b (no grid, no cut phase), tai50b's settings: the depth-3 subtree of admm's section 4, baseline arm (233 DNN nodes there)
$H/ident.sh tai30b_F3 $B -p ../experiments/bound_lab/admm/bb_tai30b_base.par -u 637117113 -F 0:3,1:7,2:10 ../qaplib/tai30b.dat
