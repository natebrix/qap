#!/bin/bash
# README.md from README.tmpl.md: each line "@@ CMD" is replaced by CMD's output (run here).
cd "$(dirname "$0")"
while IFS= read -r l; do case "$l" in "@@ "*) eval "${l#@@ }";; *) printf '%s\n' "$l";; esac; done < README.tmpl.md > README.md
