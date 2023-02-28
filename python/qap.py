# -------------------------------------------------------------------------------------------
# qap.py
# Nathan Brixius, December 2018
#
# Simple Quadratic Assignment Problem utilities.
#
# -------------------------------------------------------------------------------------------
import io
import numpy as np


def qap_obj(a, b, p):
    return sum([a[i, j] * b[p_i, p_j] for i, p_i in enumerate(p) for j, p_j in enumerate(p)])


def inv_perm(p):
    p_inv = p.copy()
    for i, p_i in enumerate(p):
        p_inv[p_i] = i
    return p_inv


def parse_perm(t, base=1, sep=' '):
    return np.array(list(map(int, t.split(sep)))) - base


def write_qaplib(name, a, b, c=None, format='%f'):
    # Writes the flow and distance matrices for a QAP in QAPLIB format
    n = len(a)
    s = io.BytesIO()
    s.write(bytes(str(n).encode("utf-8")))
    s.write(b'\n')
    s.write(b'\n')
    np.savetxt(s, a, fmt=format)
    s.write(b'\n')
    s.write(b'\n')
    np.savetxt(s, b, fmt=format)
    if c is not None:
        s.write(b'\n')
        s.write(b'\n')
        np.savetxt(s, c, fmt=format)
    t = s.getvalue().decode()
    
    with open(name, 'w') as f:
        f.write(s.getvalue().decode())
