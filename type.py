# -------------------------------------------------------------------------------------------
# type.py
# Nathan Brixius, December 2018
#
# Optimizing the arrangement of letters on a 19th century Crown typewriter. See:
#  - http://hardmath123.github.io/crown-typewriter.html
#  - https://nathanbrixius.wordpress.com/2018/11/26/optimizing-19th-century-typewriters/
#
# -------------------------------------------------------------------------------------------
import io
import re
import numpy as np

s1 = 'XQKGBPMCOFLANDTHERISUWYJVZ'
s2 = 'ZKVGWCDNIAHTESROLUMFYBPXJQ'


def get_all_docs():
    files = ['35-0.txt', 'pg1661.txt', 'pg42671.txt']
    text = ""
    regex = re.compile('[^a-zA-Z]')
    for f in files:
        t = open(f).read().lower()
        t2 = regex.sub('', t)
        text += t2
    return text


def idx(c):
    return ord(c.lower()) - ord('a')


def adj(text):
    a = np.zeros([26, 26])
    last = 0
    for i, c in enumerate(text):
        cur = idx(c)
        if i > 0:
            a[last, cur] += 1
        last = cur
    return a


def str_to_perm(s):
    return [idx(c) for c in s]


def shifts(a, b, s):
    p = str_to_perm(s)
    return qap_obj(a, b, inv_perm(p))


def dist(n):
    a = np.zeros([26, 26])
    for i in range(26):
        for j in range(26):
            a[i, j] = abs(i - j)
    return a


def get_type(p):
    ip = inv_perm(p)
    return [chr(i + ord('a')) for i in ip]

