"""Load tho40-est-dives's record (infra/phase1/runs/tho40-est-dives/): units, strata and dives, each
stratum and dive tied to its unit (units.tsv's line ranges in dives.jsonl.gz)."""
import csv, gzip, json, os

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "infra", "phase1", "runs", "tho40-est-dives")


def load(r=R):
    lines = gzip.open(os.path.join(r, "dives.jsonl.gz"), "rt").read().splitlines()
    units = list(csv.DictReader(open(os.path.join(r, "units.tsv")), delimiter="\t"))
    strata, dives, samples = {}, {}, []
    for u in units:
        for k in ("depth_min", "depth_max", "nodes", "attempt"):
            u[k] = int(u[k]) if u[k] not in ("", None) else None
        u["strata"], u["dives"] = [], 0
        if u["has_dives"] != "yes" or u["line_from"] in ("", None):
            continue
        for line in lines[int(u["line_from"]):int(u["line_to"])]:
            x = json.loads(line)
            x["unit"] = u["uid"]
            if x["type"] == "stratum":
                strata[x["key"]] = x
                u["strata"].append(x["key"])
            elif x["type"] == "dive":
                dives.setdefault(x["key"], []).append(x)
                u["dives"] += 1
            else:
                samples.append(x)
                u["sample"] = x
    for ds in dives.values():
        for x in ds:
            x["by_depth"] = {int(k): float(v) for k, v in x["by_depth"].items()}
            x["cpu_by_depth"] = {int(k): float(v) for k, v in x["cpu_by_depth"].items()}
    est = json.load(open(os.path.join(r, "estimate.json")))
    return units, strata, dives, samples, est
