"""Zenodo record 17582116 lookup helper (for future maintainers).

The committed `data/datasets.csv` already holds the exact resolved
filenames/sizes/md5s/URLs for the six k=31 eulertigs collections this
benchmark uses (cod, kestrel, human, ncbi-virus, se, hprc) plus their
matching query FASTQ -- this module is only for re-verifying those values or
extending the catalog later (e.g. to the record's k=63 files, which are out
of scope for this benchmark). Nothing in the normal `download`/`run`/`export`
path calls this module or needs network access to Zenodo.

The record's own description text does not spell out which query FASTQ
(SRR/ERR accession) belongs to which species collection -- that mapping was
resolved by hand via the ENA API (`scientific_name` for each accession) when
this benchmark's dataset catalog was built; see benchmark/sshash/README.md.
"""

from __future__ import annotations

import argparse
import json
import sys
import urllib.request

RECORD_ID = "17582116"
API_URL = f"https://zenodo.org/api/records/{RECORD_ID}"
FILE_URL_TMPL = f"https://zenodo.org/api/records/{RECORD_ID}/files/{{filename}}/content"


def fetch_record(api_url: str = API_URL, timeout: float = 30.0) -> dict:
    req = urllib.request.Request(api_url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def list_files(record: dict) -> list[dict]:
    """[{filename, size, checksum, url}, ...] from a fetch_record() result."""
    out = []
    for entry in record.get("files", []):
        filename = entry.get("key") or entry.get("filename")
        checksum = (entry.get("checksum") or "").removeprefix("md5:")
        out.append(
            {
                "filename": filename,
                "size": entry.get("size"),
                "md5": checksum,
                "url": FILE_URL_TMPL.format(filename=filename),
            }
        )
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="scripts zenodo", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list-files", help="fetch the record and print its file catalog as JSON")
    args = ap.parse_args(argv)

    if args.cmd == "list-files":
        record = fetch_record()
        files = list_files(record)
        print(json.dumps(files, indent=2))
        return 0

    print("unknown command", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
