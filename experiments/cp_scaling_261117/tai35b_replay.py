#!/usr/bin/env python3
"""Run the pooled tai35b control with the shared production-replay harness."""

import argparse
import importlib.util
from pathlib import Path
import re
from types import SimpleNamespace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("harness", type=Path)
    parser.add_argument("action", choices=["capture", "run"])
    parser.add_argument("binary", type=Path)
    parser.add_argument("--arm", default="stock")
    parser.add_argument("--cap", type=int, default=200)
    parser.add_argument("--exponent", type=float, default=0)
    parser.add_argument("--window", type=int, default=0)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location("replay", args.harness)
    replay = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(replay)

    tag = "tai35b-d2-mirror"
    replay.PURPOSES[tag] = ("tai35b", "19:12,18:26")

    def settings(instance, folder, cap=None):
        assert instance == "tai35b"
        regions = replay.REPO / "experiments/bound_lab/tai_solves/regions/tai35b.comb.txt"
        source = (replay.RECIPES / "tai50b.par").read_text()
        for key, value in {
            "dnn.depths": "0:33",
            "cut.regions": str(regions),
            "cut.cp_iters": str(args.cap if cap is None else cap),
        }.items():
            source, matches = re.subn(rf"^{re.escape(key)}\s*=[^\n]*$", f"{key} = {value}", source, flags=re.M)
            assert matches == 1, key
        recipe = folder / "params.par"
        recipe.write_text(source)
        return ["-p", str(recipe), "-u", "283315446"], replay.REPO / "qaplib/tai35b.dat", recipe

    replay.settings = settings
    if args.action == "capture":
        replay.capture(SimpleNamespace(
            tag=f"{tag}-parent", fixes="19:12", instance="tai35b", line="r18", parent=None, binary=args.binary
        ))
        replay.extract(SimpleNamespace(node=tag, source=f"{tag}-parent"))
    else:
        options = (
            [f"DNN_CP_DSCALE={args.exponent:g}", f"DNN_CP_DSCALE_UNTIL={args.window}"]
            if args.exponent else []
        )
        replay.run(SimpleNamespace(
            node=tag, arm=args.arm, binary=args.binary, cap=args.cap, env=options, complete=False
        ))


if __name__ == "__main__":
    main()
