#!/usr/bin/env python3
"""Subsample Kraken 2 library.fna files into db/standard_p10 / db/standard_p0.5.

Reads each library.fna once (optionally in parallel shards) and writes every
destination in the same pass. Keep decisions are blake2b(seed, seqid), so
results do not depend on -j or library processing order.

Usage:
  create_sample_db.py [-j JOBS] [--seed SEED] SRC_DB FRAC DST [FRAC DST ...]
"""

from __future__ import print_function

import argparse
import hashlib
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

BUF = 16 << 20
MIN_SHARD = 256 << 20
SEED_DEFAULT = 42
PROGRESS_BYTES = 512 << 20


def copy_range(out_fd, in_fd, start, size):
    if size <= 0:
        return
    remaining = size
    off = start
    if hasattr(os, "sendfile"):
        try:
            while remaining:
                n = os.sendfile(out_fd, in_fd, off, remaining)
                if n <= 0:
                    break
                off += n
                remaining -= n
            if remaining == 0:
                return
        except OSError:
            remaining = size
            off = start
    while remaining:
        chunk = os.pread(in_fd, min(remaining, BUF), off)
        if not chunk:
            break
        os.write(out_fd, chunk)
        off += len(chunk)
        remaining -= len(chunk)


def seqid_from_header(header):
    body = header[1:] if header.startswith(b">") else header
    parts = body.split()
    return parts[0] if parts else b""


def keep_score(seqid, seed):
    """Deterministic [0, 1) from seed+seqid (order- and shard-independent)."""
    h = hashlib.blake2b(digest_size=8, person=b"k2sample")
    h.update((seed & 0xFFFFFFFFFFFFFFFF).to_bytes(8, "little"))
    h.update(seqid)
    return int.from_bytes(h.digest(), "little") * (2.0**-64)


def filter_map(src_map, dst_map, kept):
    with open(src_map, "rb", buffering=BUF) as fin, open(
        dst_map, "wb", buffering=BUF
    ) as fout:
        for line in fin:
            if line.startswith(b"#") or not line.strip():
                fout.write(line)
                continue
            tab1 = line.find(b"\t")
            if tab1 < 0:
                continue
            tab2 = line.find(b"\t", tab1 + 1)
            sid = line[tab1 + 1 : tab2 if tab2 >= 0 else None]
            if sid in kept:
                fout.write(line)


def align_record_start(fd, lo, file_size):
    """File offset of the first FASTA record starting at or after lo."""
    if lo <= 0:
        return 0
    if lo >= file_size:
        return file_size
    prev = os.pread(fd, 1, lo - 1)
    cur = os.pread(fd, 1, lo)
    if prev == b"\n" and cur == b">":
        return lo
    pos = lo
    carry_nl = False
    while pos < file_size:
        chunk = os.pread(fd, BUF, pos)
        if not chunk:
            break
        if carry_nl and chunk[0] == 62:
            return pos
        idx = chunk.find(b"\n>")
        if idx >= 0:
            return pos + idx + 1
        carry_nl = chunk[-1] == 10
        pos += len(chunk)
    return file_size


def n_shards_for(size, jobs):
    if jobs <= 1 or size < MIN_SHARD:
        return 1
    return min(jobs, max(1, size // MIN_SHARD))


def sample_shard(src_fa, lo, hi, dest_paths_fracs, seed, src_size, label):
    """Copy records with start in [lo, hi) to dest files. Returns stats dict."""
    in_fd = os.open(src_fa, os.O_RDONLY)
    out_fds = []
    try:
        try:
            os.posix_fadvise(
                in_fd, max(0, lo - 1), max(0, hi - lo + BUF), os.POSIX_FADV_SEQUENTIAL
            )
        except (AttributeError, OSError, OverflowError):
            pass
        start = align_record_start(in_fd, lo, src_size)
        n_total = 0
        n_kept = [0] * len(dest_paths_fracs)
        kept_ids = [set() for _ in dest_paths_fracs]
        if start >= hi or start >= src_size:
            for path, _frac in dest_paths_fracs:
                os.close(os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644))
            return {
                "n_total": 0,
                "n_kept": n_kept,
                "kept": [[] for _ in dest_paths_fracs],
            }

        for path, _frac in dest_paths_fracs:
            out_fds.append(os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644))

        header_buf = bytearray()
        in_header = True
        rec_start = start
        last_nl = False
        keep = [False] * len(dest_paths_fracs)
        seqid = b""
        os.lseek(in_fd, start, os.SEEK_SET)
        consumed = start
        last_print = time.monotonic()
        stop = False

        def finish_header():
            nonlocal n_total, seqid, keep
            n_total += 1
            seqid = seqid_from_header(bytes(header_buf))
            score = keep_score(seqid, seed)
            keep = [score < frac for _path, frac in dest_paths_fracs]

        def flush_record(end):
            size = end - rec_start
            for i, ofd in enumerate(out_fds):
                if keep[i]:
                    copy_range(ofd, in_fd, rec_start, size)
                    kept_ids[i].add(seqid)
                    n_kept[i] += 1

        while not stop:
            chunk = os.read(in_fd, BUF)
            if not chunk:
                if rec_start < hi:
                    if not in_header:
                        flush_record(consumed)
                    elif header_buf:
                        finish_header()
                        flush_record(consumed)
                break
            chunk_off = consumed
            i = 0
            n = len(chunk)
            while i < n:
                if in_header:
                    if rec_start >= hi:
                        stop = True
                        break
                    nl = chunk.find(b"\n", i)
                    if nl < 0:
                        header_buf.extend(chunk[i:])
                        break
                    header_buf.extend(chunk[i : nl + 1])
                    finish_header()
                    in_header = False
                    i = nl + 1
                    last_nl = True
                    continue
                if last_nl and chunk[i] == 62:
                    flush_record(chunk_off + i)
                    rec_start = chunk_off + i
                    if rec_start >= hi:
                        stop = True
                        break
                    in_header = True
                    header_buf = bytearray(b">")
                    i += 1
                    last_nl = False
                    continue
                pos = chunk.find(b"\n>", i)
                if pos < 0:
                    last_nl = chunk[-1] == 10
                    break
                rec_end = chunk_off + pos + 1
                flush_record(rec_end)
                rec_start = rec_end
                if rec_start >= hi:
                    stop = True
                    break
                in_header = True
                header_buf = bytearray()
                i = pos + 1
                last_nl = False
            consumed += n
            now = time.monotonic()
            if src_size > PROGRESS_BYTES and now - last_print >= 2.0:
                last_print = now
                span = max(1, hi - lo)
                pct = 100.0 * min(consumed, hi) / src_size
                print(
                    f"    {label}  {pct:.0f}% of file  {n_total} seqs in shard",
                    flush=True,
                )
        return {
            "n_total": n_total,
            "n_kept": n_kept,
            "kept": [list(s) for s in kept_ids],
        }
    finally:
        for ofd in out_fds:
            os.close(ofd)
        os.close(in_fd)


def _shard_task(task):
    (
        name,
        shard_i,
        n_shards,
        src_fa,
        lo,
        hi,
        dest_paths_fracs,
        seed,
        src_size,
    ) = task
    label = name if n_shards == 1 else f"{name} [{shard_i + 1}/{n_shards}]"
    rec = sample_shard(src_fa, lo, hi, dest_paths_fracs, seed, src_size, label)
    rec["name"] = name
    rec["shard_i"] = shard_i
    return rec


def concat_parts(parts, dest):
    if len(parts) == 1:
        if parts[0] != dest and os.path.isfile(parts[0]):
            os.replace(parts[0], dest)
        return
    out_fd = os.open(dest, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    try:
        for p in parts:
            if not os.path.isfile(p):
                continue
            in_fd = os.open(p, os.O_RDONLY)
            try:
                copy_range(out_fd, in_fd, 0, os.fstat(in_fd).st_size)
            finally:
                os.close(in_fd)
            os.remove(p)
    finally:
        os.close(out_fd)


def list_libraries(src_db):
    src_lib = os.path.join(src_db, "library")
    libs = []
    for name in sorted(os.listdir(src_lib)):
        sdir = os.path.join(src_lib, name)
        if not os.path.isdir(sdir):
            continue
        src_fa = os.path.join(sdir, "library.fna")
        src_map = os.path.join(sdir, "prelim_map.txt")
        if not os.path.isfile(src_fa) and not os.path.isfile(src_map):
            continue
        size = os.path.getsize(src_fa) if os.path.isfile(src_fa) else 0
        libs.append((name, sdir, src_fa, src_map, size))
    return libs


def process(src_db, pairs, jobs, seed):
    jobs = max(1, int(jobs))
    libs = list_libraries(src_db)
    if not libs:
        print("  no libraries found", flush=True)
        return

    tasks = []
    meta = {}
    for name, sdir, src_fa, src_map, size in libs:
        ddirs = []
        for _frac, dst_db in pairs:
            ddir = os.path.join(dst_db, "library", name)
            os.makedirs(ddir, exist_ok=True)
            ddirs.append(ddir)
        nsh = n_shards_for(size, jobs) if os.path.isfile(src_fa) else 0
        part_paths = [[] for _ in pairs]
        if nsh:
            print(
                f"  {name}: sampling {size / 1e9:.2f} GB"
                f"{'' if nsh == 1 else f' ({nsh} shards)'}...",
                flush=True,
            )
            for i in range(nsh):
                lo = i * size // nsh
                hi = size if i + 1 == nsh else (i + 1) * size // nsh
                dest_paths_fracs = []
                for di, (frac, _dst) in enumerate(pairs):
                    path = os.path.join(
                        ddirs[di],
                        "library.fna" if nsh == 1 else f"library.fna.part{i:04d}",
                    )
                    part_paths[di].append(path)
                    dest_paths_fracs.append((path, frac))
                tasks.append(
                    (
                        name,
                        i,
                        nsh,
                        src_fa,
                        lo,
                        hi,
                        dest_paths_fracs,
                        seed,
                        size,
                    )
                )
        else:
            print(f"  {name}: no library.fna, skipping fasta", flush=True)
        meta[name] = {
            "sdir": sdir,
            "src_fa": src_fa,
            "src_map": src_map,
            "ddirs": ddirs,
            "n_shards": nsh,
            "part_paths": part_paths,
            "size": size,
        }

    results = {name: [None] * meta[name]["n_shards"] for name in meta}
    if not tasks:
        shard_recs = []
    elif jobs <= 1 or len(tasks) == 1:
        shard_recs = [_shard_task(t) for t in tasks]
    else:
        workers = min(jobs, len(tasks))
        print(f"  using {workers} workers for {len(tasks)} shard(s)", flush=True)
        shard_recs = []
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(_shard_task, t) for t in tasks]
            for fut in as_completed(futs):
                rec = fut.result()
                shard_recs.append(rec)
                m = meta[rec["name"]]
                if m["n_shards"] > 1:
                    done = sum(1 for r in shard_recs if r["name"] == rec["name"])
                    print(
                        f"    {rec['name']}: {done}/{m['n_shards']} shards",
                        flush=True,
                    )

    for rec in shard_recs:
        results[rec["name"]][rec["shard_i"]] = rec

    for name, m in meta.items():
        kept = [set() for _ in pairs]
        n_total = 0
        n_kept = [0] * len(pairs)
        if m["n_shards"]:
            for rec in results[name]:
                n_total += rec["n_total"]
                for i, ids in enumerate(rec["kept"]):
                    kept[i].update(ids)
                    n_kept[i] += rec["n_kept"][i]
            for di, (_frac, _dst) in enumerate(pairs):
                dest_fa = os.path.join(m["ddirs"][di], "library.fna")
                concat_parts(m["part_paths"][di], dest_fa)
        if len(pairs) == 1:
            print(f"  {name}: kept {n_kept[0]}/{n_total} sequences", flush=True)
        else:
            for (_frac, dst_db), nk in zip(pairs, n_kept):
                label = os.path.basename(dst_db.rstrip(os.sep))
                print(f"    {name}/{label}: kept {nk}/{n_total} sequences", flush=True)
        if os.path.isfile(m["src_map"]):
            for ddir, ids in zip(m["ddirs"], kept):
                filter_map(
                    m["src_map"], os.path.join(ddir, "prelim_map.txt"), ids
                )
        for ddir in m["ddirs"]:
            for fn in ("assembly_summary.txt", "manifest.txt"):
                sp = os.path.join(m["sdir"], fn)
                dp = os.path.join(ddir, fn)
                if os.path.isfile(sp) and not os.path.lexists(dp):
                    os.symlink(os.path.abspath(sp), dp)


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Subsample Kraken 2 library.fna files into sample DBs."
    )
    p.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 1)
    p.add_argument("--seed", type=int, default=SEED_DEFAULT)
    p.add_argument("src_db", help="Full Kraken DB (contains library/)")
    p.add_argument(
        "pairs",
        nargs="+",
        help="FRAC DST_DB pairs, e.g. 0.1 .../db/standard_p10 0.005 .../db/standard_p0.5",
    )
    args = p.parse_args(argv)
    if len(args.pairs) % 2 != 0:
        p.error("expected FRAC DST_DB pairs after SRC_DB")
    pairs = []
    for i in range(0, len(args.pairs), 2):
        pairs.append((float(args.pairs[i]), args.pairs[i + 1]))
    args.pairs = pairs
    return args


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    process(args.src_db, args.pairs, args.jobs, args.seed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
