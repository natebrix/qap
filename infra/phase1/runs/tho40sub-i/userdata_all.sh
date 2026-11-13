#!/bin/bash
export AWS_DEFAULT_REGION=eu-north-1 QAP_RUN=tho40sub-i QAP_BACKEND='aws:qap-coord-033593658551-eu-north-1?region=eu-north-1' QAP_BUCKET=qap-coord-033593658551-eu-north-1 QAP_ROLE=all
export QAP_RECIPE=tho40sub QAP_NO_UB=0 QAP_REFERENCE=1 QAP_COMMIT=dd70e71c QAP_VISIBILITY=public QAP_GATES='none'

shutdown -h +40
mkdir -p /opt/qap && cd /opt/qap
for i in 1 2 3 4 5 6; do aws s3 cp --only-show-errors s3://qap-coord-033593658551-eu-north-1/src/qap-dd70e71c.tar.gz src.tgz && break; sleep 10; done
tar xzf src.tgz && bash infra/phase1/node.sh
shutdown -h now
