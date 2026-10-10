#!/usr/bin/env python3
"""Export Wood et al. 2019 Additional file 1 tables to CSV.

Default XLSX URL (Genome Biology ESM):
  13059_2019_1891_MOESM1_ESM.xlsx

Usage:
  export_sup_tables.py --xlsx benchmark/kraken2/kraken2_sup.xlsx --out-dir benchmark/kraken2/tables
  export_sup_tables.py --download --out-dir benchmark/kraken2/tables
"""

from __future__ import print_function

import argparse
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "common"))
import units
from xlsx_util import load_xlsx, row_values
from xlsx_util import write_csv as _shared_write_csv

ESM_URL = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1186%2Fs13059-019-1891-0/MediaObjects/"
    "13059_2019_1891_MOESM1_ESM.xlsx"
)
UA = {"User-Agent": "kraken2-paper-eval"}


def _f(x):
    if x is None or x == "":
        return ""
    try:
        return float(x)
    except ValueError:
        return x


def write_csv(path, header, rows):
    _shared_write_csv(path, header, (dict(zip(header, row)) for row in rows))
    print("  wrote", path, "({} rows)".format(len(rows)))


def export_s3(sheets, out_dir):
    rows = sheets["S3 - Exclusion genomes"]
    out = []
    for r in sorted(rows):
        vals = row_values(rows[r], 1, 10)
        if not vals or not str(vals[0]).isdigit():
            continue
        taxid = int(vals[0])
        name = vals[1]
        sk, phylum, lin3, lin4, family, genus, species = vals[3:10]
        domain = "virus" if sk == "Viruses" else "prokaryote"
        out.append(
            [
                taxid,
                name,
                domain,
                sk,
                phylum,
                lin3,
                lin4,
                family,
                genus,
                species,
            ]
        )
    path = os.path.join(out_dir, "table_s3.csv")
    write_csv(
        path,
        [
            "taxid",
            "scientific_name",
            "domain",
            "superkingdom",
            "phylum",
            "lineage_3",
            "lineage_4",
            "family",
            "genus",
            "species",
        ],
        out,
    )
    n_p = sum(1 for x in out if x[2] == "prokaryote")
    n_v = sum(1 for x in out if x[2] == "virus")
    print("    {} prokaryote + {} virus = {} excluded taxa".format(n_p, n_v, len(out)))
    return path


def export_s5(sheets, out_dir):
    rows = sheets["S5 - FDA-ARGOS Evaluation"]
    out = []
    expected = []
    for r in sorted(rows):
        vals = row_values(rows[r], 1, 29)
        if not vals or vals[0] in ("", "Species", "Averages"):
            continue
        if not str(vals[1]).startswith("SRR"):
            continue
        species, run = vals[0], vals[1]
        k2_conc = _f(vals[2])
        k2_disc = _f(vals[9])
        k2_anc = _f(vals[16])
        k2_un = _f(vals[23])
        out.append(
            [
                species,
                run,
                k2_conc,
                _f(vals[3]),
                _f(vals[4]),
                _f(vals[5]),
                _f(vals[6]),
                _f(vals[7]),
                k2_disc,
                k2_anc,
                k2_un,
            ]
        )
        expected.append([species, run, k2_conc, k2_disc, k2_anc, k2_un])
    write_csv(
        os.path.join(out_dir, "table_s5.csv"),
        [
            "species",
            "sra_run",
            "k2_genus_concordant",
            "k1_genus_concordant",
            "centrifuge_genus_concordant",
            "clark_genus_concordant",
            "kaiju_genus_concordant",
            "k2x_genus_concordant",
            "k2_genus_discordant",
            "k2_genus_ancestor",
            "k2_genus_unclassified",
        ],
        out,
    )
    write_csv(
        os.path.join(out_dir, "table_s5_kraken2_expected.csv"),
        [
            "species",
            "sra_run",
            "genus_concordant",
            "genus_discordant",
            "genus_ancestor",
            "genus_unclassified",
        ],
        expected,
    )
    return os.path.join(out_dir, "table_s5.csv")


def _kb_to_gb(v):
    """Table S1 peak RSS (kbytes, KiB as GNU time reports it) -> decimal GB."""
    f = _f(v)
    return units.gb_from_kib(f) if isinstance(f, (int, float)) else f


def export_s1(sheets, out_dir):
    rows = sheets["S1 - Classifier Comparison"]
    out = []
    for r in sorted(rows):
        vals = row_values(rows[r], 1, 14)
        if not vals or vals[0] in ("", "Classifier"):
            continue
        if vals[1] not in ("nucleotide", "translated"):
            continue
        out.append(
            [
                vals[0],
                vals[1],
                _f(vals[3]),
                _f(vals[4]),
                _f(vals[5]),
                _f(vals[7]),
                _f(vals[8]),
                _f(vals[9]),
                _kb_to_gb(vals[11]),
                _f(vals[12]),
                _f(vals[13]),
            ]
        )
    write_csv(
        os.path.join(out_dir, "table_s1.csv"),
        [
            "classifier",
            "search_type",
            "prok_genus_sensitivity",
            "prok_genus_ppv",
            "prok_genus_f1",
            "viral_genus_sensitivity",
            "viral_genus_ppv",
            "viral_genus_f1",
            "max_rss_gb",
            "runtime_s",
            "mreads_per_min",
        ],
        out,
    )


def download_xlsx(dest):
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    print("  downloading", ESM_URL)
    req = urllib.request.Request(ESM_URL, headers=UA)
    with urllib.request.urlopen(req, timeout=120) as resp, open(dest, "wb") as out:
        out.write(resp.read())
    print("  saved", dest)


def parse_args(argv):
    p = argparse.ArgumentParser(description="Export Kraken 2 paper supplementary tables.")
    p.add_argument("--xlsx", default="", help="Path to Additional file 1 XLSX")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--download", action="store_true", help="Download XLSX if missing")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    xlsx = args.xlsx
    if not xlsx:
        xlsx = os.path.join(os.path.dirname(args.out_dir), "kraken2_sup.xlsx")
        if args.out_dir.rstrip("/").endswith("tables"):
            xlsx = os.path.join(os.path.dirname(args.out_dir), "kraken2_sup.xlsx")
    if not os.path.isfile(xlsx):
        if args.download:
            download_xlsx(xlsx)
        else:
            sys.exit("missing XLSX (pass --xlsx or --download): " + xlsx)
    print("[*] Reading", xlsx)
    sheets = load_xlsx(xlsx)
    os.makedirs(args.out_dir, exist_ok=True)
    export_s1(sheets, args.out_dir)
    export_s3(sheets, args.out_dir)
    export_s5(sheets, args.out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
