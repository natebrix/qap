#!/usr/bin/env python3
"""Phase 1 (README.md): the few backend operations the shell scripts need,
through coord/backend.py, so that they work alike on AWS and on a local moto
server (the local test).

  bk.py BACKEND status RUN          the run's status (running, paused, finished) or "none"
  bk.py BACKEND put KEY FILE        upload FILE to the store under KEY
  bk.py BACKEND get KEY FILE        download KEY to FILE (exit 1 if there is none)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from coord.backend import open_backend  # noqa: E402


def main(argv):
    if len(argv) < 3:
        sys.exit(__doc__)
    b = open_backend(argv[1])
    cmd = argv[2]
    if cmd == "status":
        r = b.get_run(argv[3])
        print(r.get("status", "?") if r else "none")
    elif cmd == "put":
        b.put_file(argv[3], argv[4])
    elif cmd == "get":
        return 0 if b.get_file(argv[3], argv[4]) else 1
    else:
        sys.exit(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
