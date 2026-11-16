"""Validity over every DNN trace under runs/ (tai50b, incumbent 458,821,518 = best known + 1, cutoff 458,821,517.05):
a node all of whose fixes agree with QAPSOLN/tai50b.sln contains a solution of value 458,821,517, so it must not
be fathomed and no certified value there may exceed 458,821,517.   python3 validity.py"""
import glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nodes as N
H = os.path.dirname(os.path.abspath(__file__)); BKV = 458821517
fs = [f for f in glob.glob(H + "/runs/**/dnn.txt", recursive=True) if "/t30_" not in f]
al = N.load(fs); on = [d for d in al if d["onpath"]]
print(f"{len(al)} traced node bounds in {len(fs)} traces; {len(on)} on the path of the best known solution")
bad = [d for d in on if d["fath"] or d["final"] > BKV]
seen = {}
for d in sorted(on, key=lambda d: (d["depth"], d["file"])):
    print(f"  depth {d['depth']} {','.join(d['fixes']) or '-'}: certified {d['final']:,.4f}  ({BKV - d['final']:,.0f} under)  {'FATHOMED' if d['fath'] else 'branched'}  [{os.path.relpath(d['file'], H)}]")
print("violations:", len(bad)); print("largest certified value at an unfathomed node anywhere:", f"{max(d['final'] for d in al if not d['fath']):,.4f}")
print("largest certified value at a fathomed node:", f"{max(d['final'] for d in al if d['fath']):,.4f}", "(a fathom stops at the first certificate above the cutoff)")
