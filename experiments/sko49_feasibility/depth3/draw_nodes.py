# 16 random depth-3 sko49 nodes, drawn as the wil50 sample was (NOTES.md section 2,
# depth3/run.sh): three random facilities and three random locations, paired.
import random, sys
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261102
random.seed(seed)
for k in range(16):
    f = random.sample(range(49), 3); l = random.sample(range(49), 3)
    print(",".join(f"{a}:{b}" for a, b in zip(f, l)))
