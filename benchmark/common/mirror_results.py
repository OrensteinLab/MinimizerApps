#!/usr/bin/env python3
"""Mirror a finished result set into a new one, so that later steps write there and never into the source.

    mirror_results.py SRC_RESULTS_DIR DST_RESULTS_DIR

Every file of SRC is created in DST at the same path: files of LINK_MIN_BYTES or more (databases,
indexes, read outputs: only ever read) as hard links, symbolic links as symbolic links, and all
smaller files (logs, JSON, CSVs, caches, workbooks -- the ones an export or a script may rewrite in
place) as copies, so that writing any file of DST leaves SRC unchanged. Files that already exist in DST are kept (re-running adds only what is
missing). SRC and DST must be on the same file system. Used by the profiling scripts
(benchmark/<app>/experiments/profile_run.sh). Standard library only.
"""
import os
import shutil
import sys

LINK_MIN_BYTES = 8 << 20


def mirror(src, dst):
    n_link = n_copy = 0
    for root, dirs, files in os.walk(src):
        out = os.path.join(dst, os.path.relpath(root, src))
        os.makedirs(out, exist_ok=True)
        for name in files:
            a, b = os.path.join(root, name), os.path.join(out, name)
            if os.path.lexists(b):
                continue
            if os.path.islink(a):   # a symbolic link stays one, with the same target
                os.symlink(os.readlink(a), b)
                n_copy += 1
            elif os.path.getsize(a) >= LINK_MIN_BYTES:
                os.link(a, b)
                n_link += 1
            else:
                shutil.copy2(a, b)
                n_copy += 1
    return n_link, n_copy


def main(argv):
    if len(argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    src, dst = os.path.abspath(argv[1]), os.path.abspath(argv[2])
    if not os.path.isdir(src):
        print("error: no result set %s" % src, file=sys.stderr)
        return 1
    if src == dst:
        print("error: source and destination are the same", file=sys.stderr)
        return 1
    n_link, n_copy = mirror(src, dst)
    print("[mirror] %s -> %s: %d files linked, %d copied" % (src, dst, n_link, n_copy), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
