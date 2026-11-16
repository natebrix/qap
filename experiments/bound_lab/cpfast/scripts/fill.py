#!/usr/bin/env python3
"""fill.py: replace the @@NAME@@ marks of README.md by the files fill/NAME.md (those that exist)."""
import os,re
d=os.path.dirname(os.path.abspath(__file__))+'/..'
s=open(d+'/README.md').read()
for m in sorted(set(re.findall(r'@@([A-Z0-9_]+)@@',s))):
    p=d+'/fill/'+m+'.md'
    if os.path.exists(p): s=s.replace('@@'+m+'@@',open(p).read().rstrip('\n'))
    else: print('pending',m)
open(d+'/README.md','w').write(s)
