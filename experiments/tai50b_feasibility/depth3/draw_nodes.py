# 12 random depth-3 tai50b nodes, drawn as the sko49 and wil50 samples were:
# three random facilities paired with three random locations (seed 20261108).
import random, sys
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261108
random.seed(seed)
for k in range(12):
    f = random.sample(range(50), 3); l = random.sample(range(50), 3)
    print(",".join(f"{a}:{b}" for a, b in zip(f, l)))
