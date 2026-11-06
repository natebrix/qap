#!/usr/bin/env python3
"""Compare two object files function by function: the disassembly of each
function with addresses, jump offsets and the replaced Meschach functions'
parameter lists normalized away, and each reference to a local constant
(.LCn, whose number changes when constants are added elsewhere in the file)
replaced by the constant's bytes. Prints each function that differs or is
only in one file, so a change confined to a few functions can be checked.
   python3 fndiff.py OLD.o NEW.o [-v FUNCTION_SUBSTRING]"""
import difflib
import re
import subprocess
import sys

MES = r"(m_get|m_free|m_resize|m_zero|_m_copy|v_get|v_free|v_resize|_v_copy|_in_prod)\([^)]*\)"


def run(*args):
    return subprocess.run(list(args), capture_output=True, text=True, check=True).stdout


def constants(path):
    """.LCn -> the bytes it labels (an entry of a mergeable constant section,
    or a NUL-terminated string)."""
    secs = {}
    for line in run("readelf", "-SW", path).splitlines():
        m = re.match(r"\s*\[\s*(\d+)\]\s+(\S+)\s+\S+\s+[0-9a-f]+\s+([0-9a-f]+)\s+([0-9a-f]+)\s+([0-9a-f]+)",
                     line)
        if m:
            secs[int(m.group(1))] = (m.group(2), int(m.group(3), 16), int(m.group(5), 16))
    data = open(path, "rb").read()
    lc = {}
    for line in run("readelf", "-sW", path).splitlines():
        f = line.split()
        if len(f) == 8 and f[7].startswith(".LC") and f[6].isdigit():
            name, off, es = secs[int(f[6])]
            start = off + int(f[1], 16)
            if name.startswith(".rodata.str"):
                lc[f[7]] = "str:" + data[start:data.index(b"\0", start)].decode("latin-1")
            else:
                lc[f[7]] = name + ":" + data[start:start + es].hex()
    return lc


def functions(path):
    lc = constants(path)
    out = run("objdump", "-d", "-r", "-C", "--no-show-raw-insn", path)
    funcs, name = {}, None
    for line in out.splitlines():
        m = re.match(r"^[0-9a-f]+ <(.*)>:$", line)
        if m:
            name = m.group(1)
            funcs[name] = []
            continue
        if name is None or not line.strip():
            continue
        line = re.sub(MES, r"\1", line)
        line = re.sub(r"\.LC\d+", lambda m: "[" + lc.get(m.group(0), m.group(0)) + "]", line)
        line = re.sub(r"^\s*[0-9a-f]+:\s*", "", line)                         # instruction address
        line = re.sub(r"\b[0-9a-f]+ <([^>]*?)(\+0x[0-9a-f]+)?>", r"<\1+X>", line)  # jump targets
        line = re.sub(r"#\s*<[^>]*>", "#", line)                                # rip-relative comments
        line = re.sub(r"(\t\.text[.\w]*)[+-]0x[0-9a-f]+$", r"\1+X", line)       # code moved by new code
        # a call through a relocation: the relocation on the next line names the
        # target; objdump's <label> is just whatever symbol sits at offset 0
        if line.startswith("R_X86_64") and funcs[name] and funcs[name][-1].startswith("call"):
            funcs[name][-1] = "call   <relocated>"
        funcs[name].append(line)
    return funcs


def main():
    a, b = functions(sys.argv[1]), functions(sys.argv[2])
    verbose = sys.argv[4] if len(sys.argv) > 4 and sys.argv[3] == "-v" else None
    same = 0
    for name in sorted(set(a) | set(b)):
        if name not in a:
            print(f"only in new: {name}")
        elif name not in b:
            print(f"only in old: {name}")
        elif a[name] == b[name]:
            same += 1
        else:
            print(f"differs ({len(a[name])} -> {len(b[name])} lines): {name}")
            if verbose and verbose in name:
                for d in difflib.unified_diff(a[name], b[name], lineterm="", n=2):
                    print("   ", d)
    print(f"{same} functions identical")


if __name__ == "__main__":
    main()
