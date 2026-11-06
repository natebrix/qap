#!/usr/bin/env python3
"""The profile table of README.md from logs/prof_*.log (dnntest -v -v): ms per
iteration for every piece, its share, and the achieved rate against the core's
peak where the piece is a known count of flops or bytes. Run from anywhere."""
import os, re, sys
D = os.path.dirname(os.path.abspath(__file__))
GHZ = 2.1                     # nominal clock of this Xeon (turbo unknown on the VM)
PEAK32 = 64 * GHZ             # fp32 flop/cycle on one Sapphire Rapids core (2 FMA x 16 x 2)
PEAK64 = 32 * GHZ
PEAK_AMX = 1024 * GHZ         # bf16 flop/cycle of the tile unit (16x16x32x2 per 16 cycles)

def parse(path, block=0):
    txt = open(path).read()
    if block:   # the warm child: the text after the CHILD-PASS line
        txt = txt.split("CHILD-PASS", 1)[1]
    g = lambda p: [float(x) for x in re.search(p, txt).groups()]
    prof = g(r"profile \(s\): forward ([\d.]+)  eig ([\d.]+)  expand ([\d.]+)  Y/Z-step ([\d.]+)  certify ([\d.]+)")
    eig = re.search(r"full ssyevd ([\d.]+) \((\d+) refreshes, K=(\d+)\)  filter\+RR ([\d.]+) \((\d+) steps, mean kb (\d+), mean nneg (\d+)\)", txt)
    sp = g(r"step split \(ms per subspace step\): weyl ([\d.]+)  convert ([\d.]+)  pack ([\d.]+)  W\*Y ([\d.]+)  cheb ([\d.]+)  cholqr ([\d.]+)  rr-gemm ([\d.]+)  rr-eig ([\d.]+)  x-upd ([\d.]+)")
    fr = g(r"form R \(ms per iteration\): F ([\d.]+)  syrk ([\d.]+)  symmetrize ([\d.]+)")
    res = re.search(r"RESULT .*iters=(\d+) cpu=([\d.]+)", txt)
    it = int(res.group(1)); cpu = float(res.group(2))
    r = dict(it=it, cpu=cpu, ms=1e3 * cpu / it, K=int(eig.group(3)), nref=int(eig.group(2)), refresh_s=float(eig.group(1)),
             steps=int(eig.group(5)), kb=int(eig.group(6)), nneg=int(eig.group(7)), ncert=len(re.findall(r"^it ", txt, re.M)))
    r.update(forward=1e3 * prof[0] / it, expand=1e3 * prof[2] / it, yz=1e3 * prof[3] / it, certify=1e3 * prof[4] / it)
    r.update(refresh=1e3 * r["refresh_s"] / it)
    f = r["steps"] / it   # subspace steps per iteration
    names = ["weyl", "convert", "pack", "wy", "cheb", "cholqr", "rrgemm", "rreig", "xupd"]
    for n, v in zip(names, sp): r[n] = v * f if n not in ("weyl", "convert") else v
    r.update(formF=fr[0], formsyrk=fr[1], formsym=fr[2])
    return r

rows = [("cold, sgemm (no AMX)", "prof_cold_sgemm.log", 0), ("cold, AMX", "prof_cold_amx.log", 0),
        ("cold, production (AMX, R32_NONCERT, beta schedule)", "prof_cold_prod.log", 0),
        ("warm child of 39:27,30:48 (Z + conditioned Y), production", "prof_warm_prod.log", 1),
        ("depth-4 fathom (all32 node 19, m = 46), production", "prof_d4n19_prod.log", 0)]
R = {name: parse(os.path.join(D, "logs", f), b) for name, f, b in rows}
pieces = [("W*Y products (+ pack)", lambda r: r["wy"] + r["pack"]), ("CholQR (one pass)", lambda r: r["cholqr"]),
          ("Rayleigh-Ritz: Y'T, ssyevd(kb), X = YG", lambda r: r["rrgemm"] + r["rreig"] + r["xupd"]),
          ("convert W to fp32 + Chebyshev updates", lambda r: r["convert"] + r["cheb"]),
          ("full ssyevd refresh", lambda r: r["refresh"]), ("form R (F, syrk, add-and-mirror)", lambda r: r["formF"] + r["formsyrk"] + r["formsym"]),
          ("forward V'(Y + Z/beta)V", lambda r: r["forward"]), ("expand V R V'", lambda r: r["expand"]), ("Y/Z-step", lambda r: r["yz"]),
          ("certificate (every 10)", lambda r: r["certify"])]
print("| piece | " + " | ".join(n for n, _, _ in rows) + " |")
print("|---|" + "---|" * len(rows))
for pn, fn in pieces:
    print("| %s | " % pn + " | ".join("%.1f" % fn(R[n]) for n, _, _ in rows) + " |")
print("| **total, ms per iteration** | " + " | ".join("**%.1f**" % R[n]["ms"] for n, _, _ in rows) + " |")
print("| iterations, K, mean kb / nneg | " + " | ".join("%d, %d, %d / %d" % (R[n]["it"], R[n]["K"], R[n]["kb"], R[n]["nneg"]) for n, _, _ in rows) + " |")
print()
# achieved rates on the production cold run
r = R[rows[2][0]]; K = r["K"]; kb = r["kb"]; nneg = r["nneg"]; f = r["steps"] / r["it"]
prod_ms = r["wy"] / f / 4          # one W*Y product (degree 3: 3 filter products + 1 Rayleigh-Ritz)
gf = 2.0 * K * K * kb / 1e9
print("achieved rates, production cold run (ms are per call; peaks at the nominal %.1f GHz: fp32 %.0f GF/s, fp64 %.0f GF/s, AMX bf16 %.0f GF/s):" % (GHZ, PEAK32, PEAK64, PEAK_AMX))
print("- one W*Y product on AMX: %.1f ms for %.1f GF fp32-equivalent = %.0f GF/s-equivalent = %.0f GF/s bf16 over the three tile products (%.0f%% of the tile peak)" % (prod_ms, gf, gf / prod_ms * 1e3, 3 * gf / prod_ms * 1e3, 300 * gf / prod_ms * 1e3 / PEAK_AMX))
rs = R[rows[0][0]]; prod32 = rs["wy"] / (rs["steps"] / rs["it"]) / 4
print("- one W*Y product with OpenBLAS sgemm: %.1f ms = %.0f GF/s (%.0f%% of the fp32 peak at %.1f GHz)" % (prod32, gf / prod32 * 1e3, 100 * gf / prod32 * 1e3 / PEAK32, GHZ))
chol = r["cholqr"] / f; gfc = 2.0 * K * kb * kb / 1e9
print("- CholQR (ssyrk + spotrf + strsm, %.1f GF): %.1f ms = %.0f GF/s" % (gfc, chol, gfc / chol * 1e3))
rg = r["rrgemm"] / f; gfr = 2.0 * K * kb * kb / 1e9
print("- Rayleigh-Ritz Y'T and X = YG (sgemm, %.1f GF each): %.1f and %.1f ms = %.0f GF/s" % (gfr, rg, r["xupd"] / f, gfr / rg * 1e3))
re_ = r["rreig"] / f; gfe = 10.0 / 3 * kb ** 3 / 1e9
print("- ssyevd(kb = %d) (about %.2f GF): %.1f ms = %.0f GF/s" % (kb, gfe, re_, gfe / re_ * 1e3))
ref = r["refresh_s"] / r["nref"] * 1e3; gff = (4.0 / 3 + 2 + 1) * K ** 3 / 1e9
print("- full ssyevd(K = %d) (about %.0f GF): %.0f ms per refresh = %.0f GF/s; %d refreshes in %d iterations" % (K, gff, ref, gff / ref * 1e3, r["nref"], r["it"]))
gfs = 1.0 * K * K * nneg / 1e9
print("- form R rank term, fp32 ssyrk (%.1f GF): %.1f ms = %.0f GF/s; fp64 dsyrk (cold AMX run): %.1f ms = %.0f GF/s (fp64 peak %.0f)" % (gfs, r["formsyrk"], gfs / r["formsyrk"] * 1e3, R[rows[1][0]]["formsyrk"], gfs / R[rows[1][0]]["formsyrk"] * 1e3, PEAK64))
n1 = int(K ** 0.5 - 1 + 1) ** 2 + 1
m = int(round((K - 1) ** 0.5)) + 1; n1 = m * m + 1
MB = lambda words8: words8 * 8 / 1e6
fwd_mb = MB(2 * n1 * n1) + MB(n1 * n1) + MB(2 * n1 * n1) + MB(K * K) + MB(K * K) + K * K * 4 / 1e6 + MB(2 * K * K)
print("- forward: about %.0f MB of traffic (Y, Z read; T written and rewritten by rows_D; the gather into W and Wf; the symmetrisation) in %.1f ms = %.1f GB/s" % (fwd_mb, r["forward"], fwd_mb / r["forward"]))
exp_mb = MB(n1 * n1) + MB(2 * n1 * n1) + MB(2 * n1 * n1)
print("- expand: about %.0f MB (write, rows_D, symmetrisation) in %.1f ms = %.1f GB/s" % (exp_mb, r["expand"], exp_mb / r["expand"]))
yz_mb = MB(6 * n1 * n1) + MB(n1 * n1) / 10
print("- Y/Z-step: about %.0f MB (L, VRV, Y, Z read; Y, Z written; Zc on certificate iterations) in %.1f ms = %.1f GB/s" % (yz_mb, r["yz"], yz_mb / r["yz"]))
print("- certificate: %.0f ms each (%d in %d iterations): forward %.0f + dpotrf(K) %.1f GF fp64 + Gilmore-Lawler LAPs" % (r["certify"] * r["it"] / r["ncert"], r["ncert"], r["it"], r["forward"], K ** 3 / 3 / 1e9))
