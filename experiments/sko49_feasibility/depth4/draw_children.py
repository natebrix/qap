# depth-4 children of two depth-3 misses, the way the wil50 depth-4 sample was drawn
# (../../wil50_feasibility/depth4/run.sh): one random free facility per parent, K random
# free locations. usage: draw_children.py "name=spec" "name=spec" [K] [seed]
import random, sys
args = [a for a in sys.argv[1:] if "=" in a]
nums = [int(a) for a in sys.argv[1:] if "=" not in a]
K = nums[0] if nums else 8
random.seed(nums[1] if len(nums) > 1 else 20261102)
for a in args:
    name, spec = a.split("=")
    used_f = {int(x.split(":")[0]) for x in spec.split(",")}
    used_l = {int(x.split(":")[1]) for x in spec.split(",")}
    f = random.choice([x for x in range(49) if x not in used_f])
    for l in random.sample([x for x in range(49) if x not in used_l], K):
        print(f"{name} {spec},{f}:{l}")
