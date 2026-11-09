#!/usr/bin/env python3
"""Old environment recipes into a params file's [settings] section
(docs/settings_plan.md, Phase 2).

Before Phase 2 a DNN configuration was a list of environment variables
(QAP_DNN_BETA=0.3 QAP_DNN_CUT_ZONE=0.005 ...). The solver now reads those
settings from keys in the params file (dnn.beta, cut.zone, ...) and rejects
the old names. This module turns such a list into a params file: the base
file's blocks, then a [settings] section with the keys.

    python3 bench/envpar.py latest/PARAM/nug_deep_xbar9_all.par -o run.par \\
        --env QAP_DNN_BETA=0.3 --env QAP_DNN_CUT_ZONE=0.005 --env DNN_AMX=1

    # a bash launcher that still exports the old names: write the file from the
    # process environment, then drop the converted names from it
    eval "$(python3 bench/envpar.py base.par -o run.par --from-env --shell)"

Values can also be given by key (--env cut.zone=0.005), and a base file that
has a [settings] section already is merged (the new values win). Names the map
does not know stay in the environment: the experiment switches (DNN_R32,
DNN_ADMM_WARM_BETA, ...), the diagnostics (QAP_DNN_TRACE, QAP_DNN_POOL_CHECK),
the thread variables; --shell leaves them exported, write_par() returns them.
The map is `bnbdrv --settings-map` (latest/settings.C), so it is never stale.
"""
import argparse
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOLVER = os.path.join(ROOT, "latest", "bnbdrv")
SOLVER_PREFIXES = ("QAP_", "DNN_")

_map_cache = {}


def settings_map(solver=SOLVER):
    """{old environment name: key} from the solver"""
    if solver not in _map_cache:
        out = subprocess.run([solver, "--settings-map"], capture_output=True, text=True, check=True).stdout
        _map_cache[solver] = dict(line.split() for line in out.splitlines() if len(line.split()) == 2)
    return _map_cache[solver]


def split_env(env, smap):
    """(settings, residual): the entries of env that are settings keys (by key or
    by old name) in the order given, and the rest"""
    settings, residual = {}, {}
    for k, v in env.items():
        if k in smap:
            settings[smap[k]] = v
        elif "." in k and not k.startswith(SOLVER_PREFIXES):
            settings[k] = v
        else:
            residual[k] = v
    return settings, residual


def read_par(path):
    """(body, settings): the file's text before its [settings] section and the
    section's key -> value (an empty dict when there is none)"""
    body, settings, in_settings = [], {}, False
    for line in open(path):
        s = line.strip()
        if s.lower() == "[settings]":
            in_settings = True
            continue
        if in_settings:
            if s.startswith("[") and s.endswith("]"):
                raise SystemExit(f"{path}: a {s} block after [settings] is not supported")
            s = s.split("#", 1)[0].strip()
            if s:
                k, _, v = s.partition("=")
                settings[k.strip()] = v.strip()
        else:
            body.append(line)
    return "".join(body), settings


def par_text(body, settings, comment=None):
    out = body.rstrip("\n") + "\n\n[settings]\n"
    if comment:
        out += "".join(f"# {c}\n" for c in comment.splitlines())
    out += "".join(f"{k} = {v}\n" for k, v in settings.items())
    return out


def write_par(base, env, out, solver=SOLVER, comment=None):
    """the base file plus env's settings as [settings] keys, written to out;
    returns the entries of env that are not settings (keep them in the
    environment)"""
    settings, residual = split_env(env, settings_map(solver))
    body, old = read_par(base)
    old.update(settings)
    with open(out, "w") as f:
        f.write(par_text(body, old, comment))
    return residual


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], epilog="see the module docstring")
    ap.add_argument("base", help="the params file whose blocks the output keeps")
    ap.add_argument("-o", "--out", help="the params file to write (default: stdout)")
    ap.add_argument("--env", action="append", default=[], help="NAME=VALUE, an old name or a key (repeat)")
    ap.add_argument("--from-env", action="store_true", help="also every QAP_*/DNN_* variable of the environment")
    ap.add_argument("--shell", action="store_true",
                    help="print 'unset NAME ...' for the converted environment names, for eval (needs -o)")
    ap.add_argument("--solver", default=SOLVER, help="the bnbdrv whose --settings-map to use")
    ap.add_argument("--comment", help="a comment line at the top of the section")
    a = ap.parse_args()
    env = {}
    if a.from_env:
        env.update({k: v for k, v in os.environ.items() if k.startswith(SOLVER_PREFIXES)})
    for kv in a.env:
        k, _, v = kv.partition("=")
        env[k] = v
    smap = settings_map(a.solver)
    settings, residual = split_env(env, smap)
    body, old = read_par(a.base)
    old.update(settings)
    text = par_text(body, old, a.comment)
    if a.out:
        with open(a.out, "w") as f:
            f.write(text)
    else:
        sys.stdout.write(text)
    if a.shell:
        if not a.out:
            sys.exit("--shell needs -o")
        converted = sorted(k for k in env if k in smap or ("." in k and not k.startswith(SOLVER_PREFIXES)))
        if converted:
            print("unset " + " ".join(k for k in converted if "." not in k))
    kept = sorted(k for k in residual if k.startswith(SOLVER_PREFIXES))
    if kept:
        print("envpar: kept in the environment (not settings keys): " + " ".join(kept), file=sys.stderr)


if __name__ == "__main__":
    main()
