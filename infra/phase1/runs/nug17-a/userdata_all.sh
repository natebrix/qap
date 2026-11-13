#!/bin/bash
export AWS_DEFAULT_REGION=eu-north-1 QAP_RUN=nug17-a QAP_BACKEND='aws:qap-coord-033593658551-eu-north-1?region=eu-north-1' QAP_BUCKET=qap-coord-033593658551-eu-north-1 QAP_ROLE=all
export QAP_RECIPE=nug17 QAP_NO_UB=0 QAP_REFERENCE=1 QAP_COMMIT=f12bac23 QAP_VISIBILITY=public
shutdown -h +120
mkdir -p /opt/qap && cd /opt/qap
for i in 1 2 3 4 5 6; do aws s3 cp --only-show-errors s3://qap-coord-033593658551-eu-north-1/src/qap-f12bac23.tar.gz src.tgz && break; sleep 10; done
tar xzf src.tgz && bash infra/phase1/node.sh
shutdown -h now
