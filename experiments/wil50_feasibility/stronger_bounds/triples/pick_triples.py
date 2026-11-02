#!/usr/bin/env python3
"""Random triples of local facility indices (0..m-1) for triplet_separate.py.
  python3 pick_triples.py m count [seed]   -> one 'a,b,c' per line"""
import random, sys
m, k = int(sys.argv[1]), int(sys.argv[2]); rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
seen = set()
while len(seen) < k:
    t = tuple(sorted(rng.sample(range(m), 3)))
    if t not in seen: seen.add(t); print(",".join(map(str, t)))
