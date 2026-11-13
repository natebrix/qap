#!/bin/bash
export AWS_DEFAULT_REGION=eu-north-1 QAP_RUN=nug30-a QAP_BACKEND='aws:qap-coord-033593658551-eu-north-1?region=eu-north-1' QAP_BUCKET=qap-coord-033593658551-eu-north-1 QAP_ROLE=all
export QAP_RECIPE=nug30 QAP_NO_UB=0 QAP_REFERENCE=1 QAP_COMMIT=e0e1bc2c QAP_VISIBILITY=public

shutdown -h +90
mkdir -p /opt/qap && cd /opt/qap
for i in 1 2 3 4 5 6; do aws s3 cp --only-show-errors s3://qap-coord-033593658551-eu-north-1/src/qap-e0e1bc2c.tar.gz src.tgz && break; sleep 10; done
tar xzf src.tgz && bash infra/phase1/node.sh
shutdown -h now
