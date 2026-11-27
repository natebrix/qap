# Phase 0, 20261115-2317-eu-north-1

| instance | CPU | arch | cores | AMX | ms/iteration alone | loaded | loaded / alone | CPU s per node (loaded) | market, $/h | $/core-hour | $ per 1,000 DNN nodes | gates | fathomed (reference) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| c7a.4xlarge | AMD EPYC 9R14 | x86_64 | 16 | 0 | 337.1 | 385.0 | 1.14 | 106.6 | spot 0.1225 | 0.0077 | 0.227 | all pass | 18 (18), 11 differ |
| c7i.4xlarge | Intel(R) Xeon(R) Platinum 8488C | x86_64 | 8 | 1 | 195.0 | 212.3 | 1.09 | 58.8 | spot 0.1258 | 0.0157 | 0.257 | all pass | 18 (18), 11 differ |

Per-iteration speed under load, relative to the fastest: c7a.4xlarge 1.81x, c7i.4xlarge 1.00x
- c7a.4xlarge: node 11: fathomed 290 -> fathomed 284
- c7a.4xlarge: node 17: fathomed 210 -> fathomed 205
- c7a.4xlarge: node 18: fathomed 270 -> fathomed 269
- c7a.4xlarge: node 19: fathomed 130 -> fathomed 122
- c7a.4xlarge: node 20: fathomed 170 -> fathomed 162
- c7a.4xlarge: node 21: fathomed 180 -> fathomed 178
- c7a.4xlarge: node 22: fathomed 400 -> fathomed 391
- c7a.4xlarge: node 25: fathomed 150 -> fathomed 142
- c7a.4xlarge: node 27: fathomed 350 -> fathomed 347
- c7a.4xlarge: node 31: fathomed 240 -> fathomed 233
- c7a.4xlarge: node 32: fathomed 230 -> fathomed 225
- c7i.4xlarge: node 11: fathomed 290 -> fathomed 284
- c7i.4xlarge: node 18: fathomed 270 -> fathomed 269
- c7i.4xlarge: node 19: fathomed 130 -> fathomed 122
- c7i.4xlarge: node 20: fathomed 170 -> fathomed 162
- c7i.4xlarge: node 21: fathomed 180 -> fathomed 177
- c7i.4xlarge: node 22: fathomed 400 -> fathomed 391
- c7i.4xlarge: node 25: fathomed 150 -> fathomed 142
- c7i.4xlarge: node 27: fathomed 350 -> fathomed 347
- c7i.4xlarge: node 28: fathomed 160 -> fathomed 159
- c7i.4xlarge: node 31: fathomed 240 -> fathomed 233
- c7i.4xlarge: node 32: fathomed 230 -> fathomed 225
