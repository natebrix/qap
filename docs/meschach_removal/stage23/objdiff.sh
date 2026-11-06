#!/bin/sh
# Compare the machine code of every object in two builds of latest/: the
# disassembly (objdump -d -r, demangled) with the replaced Meschach functions'
# C++ parameter lists stripped, so a call to Meschach's m_get and a call to
# mat.C's m_get(int, int) compare equal. Prints, per object, "same" or the
# number of differing lines.
#   sh objdiff.sh OLD_LATEST_DIR NEW_LATEST_DIR
old=$1 new=$2
tmp=$(mktemp -d)
norm() {
  objdump -d -r -C --no-show-raw-insn "$1" | tail -n +4 |
    sed -E 's/\b(m_get|m_free|m_resize|m_zero|_m_copy|v_get|v_free|v_resize|_v_copy|_in_prod)\([^)]*\)/\1/g'
}
for o in "$old"/*.o "$new"/*.o; do basename "$o"; done | sort -u > "$tmp/list"
while read -r b; do
  if [ ! -f "$old/$b" ]; then echo "$b: only in new"; continue; fi
  if [ ! -f "$new/$b" ]; then echo "$b: only in old"; continue; fi
  norm "$old/$b" > "$tmp/a"; norm "$new/$b" > "$tmp/b"
  if cmp -s "$tmp/a" "$tmp/b"; then echo "$b: same"
  else echo "$b: $(diff "$tmp/a" "$tmp/b" | grep -c '^[<>]') lines differ"; fi
done < "$tmp/list"
rm -rf "$tmp"
