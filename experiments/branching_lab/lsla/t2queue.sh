#!/bin/bash
# the tai35b arms of T2, one after another (2 slots each)
D=$(cd "$(dirname "$0")"; pwd)
P=${1:-tai35b}
$D/t2.sh $P prod
$D/t2.sh $P a QAP_LS=a:2:0:2
$D/t2.sh $P b2 QAP_LS=b:2:0:2 QAP_LS_CAP=30
$D/t2.sh $P b4 QAP_LS=b:4:0:2 QAP_LS_CAP=30
