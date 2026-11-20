#!/bin/bash
cd "$(dirname "$0")"
for par in - 19:26,18:12 19:12,18:26 19:26,18:8 19:26; do
  INST=tai35b THREADS=2 python3 t2.py $par 1e-4 2>&1 | grep -v "^Running\|^Includes"
done
for par in 19:26,18:12 -; do
  INST=tai35b python3 t2.py $par 0 ipm 2>&1 | grep -v "^Running\|^Includes"
done
