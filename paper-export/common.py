"""Shared definitions of the paper scripts: the result sets, the schemes and datasets the paper reports, and their names.

A result set is the folder a benchmark writes (benchmark/<method>/README.md); the paper scripts read its exported CSVs,
<results>/csvs/*.csv, and compute the statistics with benchmark/common (significance.py, density_effect.py).
"""

import csv
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "benchmark" / "common"))
from density_effect import num as fl  # noqa: E402,F401  (a CSV value as a float, None when empty or not a number)

APPS = ("kraken2", "minimap2", "kmc", "sshash")
METHOD = {"kraken2": "Kraken2", "minimap2": "Minimap2", "kmc": "KMC3", "sshash": "SSHash"}

# The polarity of Phi of every GreedyMini-based scheme the paper reports, chosen arbitrarily: A,T -> 0, C,G -> 1.
POLARITY = "at0-cg1"
POLARITIES = ("at0-cg1", "at1-cg0")

# The schemes of the paper, in the order of the tables and figures, as (scheme, order); the default order is named per method.
LABELS = ("default", "xor", "xor-pcg-hash", "GreedyMini", "mod-minimizer (xor)", "mod-minimizer (GreedyMini)")
SCHEME_ORDER = {"default": ("minimizer", None), "xor": ("minimizer", "xor"), "xor-pcg-hash": ("minimizer", "xor-pcg-hash"),
                "GreedyMini": ("minimizer", "GreedyMini"), "mod-minimizer (xor)": ("mod-minimizer", "xor"),
                "mod-minimizer (GreedyMini)": ("mod-minimizer", "GreedyMini")}
DEFAULT_ORDER = {"kraken2": "xor (default)", "minimap2": "invertible integer hash (default)", "kmc": "signature (default)",
                 "sshash": "mixer_64 (default)"}
GM_SCHEMES = {"greedymini": "GreedyMini", "mod-minimizer-greedymini": "mod-minimizer (GreedyMini)"}

# The datasets of each method, in the paper's order: (dataset id in the result set, name).
DATASETS = {
    "kraken2": [("strain-exclusion", "strain-exclusion"), ("cami", "CAMI II"), ("fda-argos", "FDA-ARGOS")],
    "minimap2": [("chm13_chrX", "CHM13 chrX"), ("grch38_wg", "GRCh38"), ("chm13_rel3_ont", "ONT reads on GRCh38")],
    "kmc": [("f_vesca", "F. vesca"), ("g_gallus", "G. gallus"), ("m_balbisiana", "M. balbisiana"),
            ("h_sapiens_2", "H. sapiens 2"), ("h_sapiens_3", "H. sapiens 3")],
    "sshash": [("cod", "cod"), ("kestrel", "kestrel"), ("human", "human"), ("ncbi-virus", "NCBI viruses"),
               ("se", "S. enterica"), ("hprc", "HPRC")],
}
KRAKEN2_LIBRARY = "strain-exclusion library"   # Kraken2 builds one database, from this library

# (k, w) of each method; Minimap2 uses its presets' k: 19 for the CLR reads, 15 for the ONT reads.
KW = {"kraken2": (31, 5), "minimap2": (19, 10), "kmc": (9, 20), "sshash": (21, 11)}
KW_DATASET = {("minimap2", "chm13_rel3_ont"): (15, 10)}


def kw(app, dataset=None):
    return KW_DATASET.get((app, dataset), KW[app])


def dataset_name(app, dataset):
    return dict(DATASETS[app]).get(dataset, dataset)


class ResultSets:
    """The result set of each method: APP=DIR from the command line."""

    def __init__(self, items):
        self.dirs = {}
        for item in items:
            app, _, d = item.partition("=")
            if app not in APPS or not d:
                raise ValueError("expected APP=DIR with APP in %s, got %r" % (", ".join(APPS), item))
            self.dirs[app] = Path(d)
        missing = [a for a in APPS if a not in self.dirs]
        if missing:
            raise ValueError("missing result sets: " + ", ".join(missing))

    def __getitem__(self, app):
        return self.dirs[app]

    def csv(self, app, name):
        """The rows of <results>/csvs/<name>."""
        with open(self.dirs[app] / "csvs" / name, newline="") as f:
            return list(csv.DictReader(f))


def label(row, app):
    """The paper's label of a result row (LABELS), or None if the paper does not report the row: GreedyMini-based schemes
    only with POLARITY; SSHash's GreedyMini only with the order on 15-mers; KMC3's only with the binary table."""
    s = row.get("scheme", "")
    if s in ("native", "default") or (app == "kraken2" and s == "xor"):
        return "default"
    if s in ("xor", "xor-pcg-hash"):
        return s
    if s == "mod-minimizer-xor":
        return "mod-minimizer (xor)"
    if s in GM_SCHEMES and row.get("projection", "") != POLARITY:
        return None
    if s == "greedymini":
        if app == "sshash" and row.get("gm_table") != "w11_k15_canonical":
            return None
        if app == "kmc" and row.get("gm_lut") == "dna":
            return None
        return "GreedyMini"
    if s == "mod-minimizer-greedymini":
        return "mod-minimizer (GreedyMini)"
    return None


def label_key(lab):
    return LABELS.index(lab)


def scheme_order(lab, app):
    """The (scheme, order) columns of a label."""
    scheme, order = SCHEME_ORDER[lab]
    return scheme, order or DEFAULT_ORDER[app]


def key_label(key, app):
    """The label of a scheme key of benchmark/common (scheme/gm_k/projection/gm_table/gm_lut/mod_t, empty parts left out,
    e.g. 'greedymini/15/at0-cg1/w11_k15_canonical')."""
    f = key.split("/")
    row = {"scheme": f[0], "projection": next((p for p in f if p.startswith("at")), ""),
           "gm_table": next((p for p in f if p.startswith("w")), ""), "gm_lut": next((p for p in f if p in ("dna", "binary")), "")}
    return label(row, app)


def minimap2_blocks(rows):
    """The rows of a Minimap2 CSV (runs or summary, in file order) as (fixed, derived): the block with the occurrence
    threshold fixed to that of the default scheme, and the block with the threshold derived from each index. The blocks
    are told apart by position (benchmark/common/scheme_tests.minimap2_fixed_block), since a derived threshold may equal
    the default scheme's."""
    from scheme_tests import minimap2_fixed_block
    in_fixed, _ = minimap2_fixed_block(rows)
    return [r for r, f in zip(rows, in_fixed) if f], [r for r, f in zip(rows, in_fixed) if not f]


def default_thresholds(R):
    """Minimap2's occurrence threshold of the default scheme on each dataset, to which the paper fixes the threshold of
    every scheme."""
    return {r["dataset"]: r["mid_occ_threshold"] for r in R.csv("minimap2", "minimap2_summary.csv") if r["scheme"] == "native"}
