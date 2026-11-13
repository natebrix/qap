# pairs.py RUNS: C against Cw on the eleven anchor nodes beyond 234 short: totals, ratio, difference in SEs,
# and the draws at levels 4 and 5 that fathom at their first solve. RUNS holds the probe.py runs (probes,
# probes2, probes3, zoneA-zoneD), whose L*/draws.json and t*/result.json are not committed.
import json, os, sys, glob, re
S = sys.argv[1]
pairs = [("Cp232","probes","Zp232","zoneD"),("Cp268","probes","Zp268","zoneD"),("Cr283","probes3","Zr283","zoneD"),
         ("Cp303","probes","Zp303","zoneD"),("Cr311","probes2","Zr311","zoneD"),("Cp326","probes","Zp326","zoneB"),
         ("Cr338","probes2","Zr338","zoneD"),("Cr350","probes3","Zr350","zoneD"),("Cr368","probes2","Zr368","zoneC"),
         ("Cr385","probes2","Zr385","zoneD"),("Cp396","probes","Zp396","zoneA")]
def tot(d):
    t = open(os.path.join(d, "summary.txt")).readline()
    m = re.search(r"DNN CPU ([\d,]+) s \(SE ([\d,]+) s", t)
    return float(m.group(1).replace(",", "")), float(m.group(2).replace(",", ""))
def lvl(d, L):
    p = os.path.join(d, f"L{L}", "draws.json")
    if not os.path.exists(p): return None
    dr = json.load(open(p)); n = len(dr["draws"]); ff = ph = fa = 0
    for j in range(n):
        doc = json.load(open(os.path.join(d, f"L{L}", f"t{j}", "result.json"))); tr = doc["trace"]
        if tr and tr[0]["cutoff"] - tr[0]["lb"] <= 0: ff += 1
        ph += int(any(t["cut"] for t in tr)); fa += int(not doc["children"])
    return n, ff, ph, fa
import math
lr = []
for c, cd, z, zd in pairs:
    dc, dz = glob.glob(f"{S}/{cd}/{c}")[0], glob.glob(f"{S}/{zd}/{z}")[0]
    (a, sa), (b, sb) = tot(dc), tot(dz)
    r = b / a; lr.append(math.log(r))
    se = (b - a) / math.sqrt(sa**2 + sb**2)
    l4c, l4z, l5c, l5z = lvl(dc, 4), lvl(dz, 4), lvl(dc, 5), lvl(dz, 5)
    print(f"{c:6} {a:10,.0f} {z:6} {b:9,.0f} ratio {r:5.2f} diff/SE {se:+5.1f} | L4 ff C {l4c[1]}/{l4c[0]} Cw {l4z[1]}/{l4z[0]} | L5 ff C {l5c and l5c[1]}/{l5c and l5c[0]} Cw {l5z and l5z[1]}/{l5z and l5z[0]}")
print("geomean ratio %.2f" % math.exp(sum(lr) / len(lr)))
