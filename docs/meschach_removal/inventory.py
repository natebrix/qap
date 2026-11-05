#!/usr/bin/env python3
"""Count uses of Meschach types, functions and macros outside src/MESCHACH.

Run from the repo root:  python3 docs/meschach_removal/inventory.py
Prints one table per category (identifier -> per-file counts) and a per-file
total. Comments are not stripped, so a commented-out call counts too; the
report says where that matters.
"""
import os
import re
import sys
from collections import defaultdict

ROOT = os.getcwd()

CATEGORIES = {
    "types": ["MAT", "VEC", "PERM", "IVEC", "BAND", "ZMAT", "ZVEC", "SPMAT"],
    "alloc/free": ["m_get", "v_get", "px_get", "iv_get", "m_free", "v_free",
                   "px_free", "iv_free", "M_FREE", "V_FREE", "PX_FREE",
                   "IV_FREE", "m_resize", "v_resize", "px_resize",
                   "iv_resize", "mem_stat_reg", "mem_stat_reg_list",
                   "MEM_STAT_REG", "mem_stat_mark", "mem_stat_free",
                   "m_get_vars", "v_get_vars"],
    "copy/init": ["m_copy", "v_copy", "px_copy", "iv_copy", "_m_copy",
                  "_v_copy", "m_zero", "v_zero", "m_ident", "px_ident",
                  "m_rand", "v_rand", "m_ones", "v_ones", "set_col",
                  "set_row", "get_col", "get_row", "m_move", "v_move",
                  "sub_mat", "sub_vec"],
    "arith (BLAS-1/2/3-like)": ["in_prod", "_in_prod", "v_add", "v_sub",
                                "sv_mlt", "v_mltadd", "m_add", "m_sub",
                                "sm_mlt", "m_mlt", "mmtr_mlt", "mtrm_mlt",
                                "mv_mlt", "vm_mlt", "m_transp", "v_norm2",
                                "v_norm_inf", "_v_norm_inf", "v_norm1",
                                "m_norm_inf", "m_norm1", "m_norm_frob",
                                "v_star", "v_sum", "v_sort", "px_rows",
                                "px_cols", "px_vec", "pxinv_vec", "ms_mltadd",
                                "v_lincomb", "v_linlist"],
    "decompositions/eigen": ["Hfactor", "makeHQ", "makeH", "trieig",
                             "symmeig", "schur", "svd", "bisvd", "LUfactor",
                             "LUsolve", "CHfactor", "CHsolve", "QRfactor",
                             "QRsolve", "BKPfactor", "Usolve", "Lsolve",
                             "hhvec", "hhtrcols", "hhtrrows"],
    "I/O and errors": ["m_output", "v_output", "px_output", "m_foutput",
                       "v_foutput", "px_foutput", "m_finput", "v_finput",
                       "m_input", "v_input", "error", "warning", "ev_err",
                       "set_err_flag", "catch", "tracecatch", "E_SIZES",
                       "E_NULL", "E_MEM", "E_SQUARE", "E_RANGE",
                       "E_INTERN", "MACHEPS"],
}

FIELD_RE = {
    "->me / .me (MAT rows)": re.compile(r"(->|\.)me\b"),
    "->ve / .ve (VEC data)": re.compile(r"(->|\.)ve\b"),
    "->pe / .pe (PERM data)": re.compile(r"(->|\.)pe\b"),
    "->base (MAT storage)": re.compile(r"->base\b"),
    "->max_m/max_n/max_size": re.compile(r"->max_(m|n|size|dim)\b"),
}

SRC_EXT = (".C", ".c", ".h", ".cpp", ".cc")


def source_files():
    for d, subs, files in os.walk(ROOT):
        rel = os.path.relpath(d, ROOT)
        # docs/ holds this study's own prototypes
        if rel.startswith(os.path.join("src", "MESCHACH")) or "/.worktrees" in d \
                or rel.startswith(".git") or rel.startswith("bench/.worktrees") \
                or rel.startswith("docs"):
            subs[:] = []
            continue
        for f in files:
            if f.endswith(SRC_EXT):
                yield os.path.relpath(os.path.join(d, f), ROOT)


def main():
    files = sorted(source_files())
    text = {f: open(os.path.join(ROOT, f), errors="replace").read() for f in files}
    per_file_total = defaultdict(int)
    for cat, names in CATEGORIES.items():
        rows = []
        for name in names:
            # "error" and "catch" are too generic: count only Meschach call
            # shapes error(E_..., ...) / catch(E_...
            if name in ("error", "warning"):
                pat = re.compile(r"\b%s\s*\(\s*(E_|WARN_)" % name)
            elif name == "catch":
                pat = re.compile(r"\bcatch\s*\(\s*E_")
            else:
                pat = re.compile(r"\b%s\b" % re.escape(name))
            counts = {f: len(pat.findall(t)) for f, t in text.items()}
            counts = {f: c for f, c in counts.items() if c}
            if counts:
                rows.append((name, counts))
                for f, c in counts.items():
                    per_file_total[f] += c
        print("## %s" % cat)
        for name, counts in rows:
            tot = sum(counts.values())
            detail = ", ".join("%s %d" % (f, c) for f, c in sorted(counts.items(), key=lambda x: -x[1]))
            print("  %-22s %5d   %s" % (name, tot, detail))
        print()
    print("## field accesses")
    for label, pat in FIELD_RE.items():
        counts = {f: len(pat.findall(t)) for f, t in text.items()}
        counts = {f: c for f, c in counts.items() if c}
        tot = sum(counts.values())
        detail = ", ".join("%s %d" % (f, c) for f, c in sorted(counts.items(), key=lambda x: -x[1]))
        print("  %-26s %5d   %s" % (label, tot, detail))
        for f, c in counts.items():
            per_file_total[f] += c
    print()
    print("## total Meschach identifiers per file")
    for f, c in sorted(per_file_total.items(), key=lambda x: -x[1]):
        print("  %5d  %s" % (c, f))


if __name__ == "__main__":
    sys.exit(main())
