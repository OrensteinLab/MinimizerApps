"""Shared layout of the `orders` sheet (and <app>_orders.csv) of all four apps.

One row per (order variant, k, w, L):

  order, k, w, L, gm_k, projection, gm_table, [gm_lut], mod_t,   <- identity
  description, cli,                                              <- what it is
  <density columns>                                              <- app-specific

with the unified parameter names

  k  minimizer length   (kmc: p, kraken2: l, sshash: m, minimap2: k)
  w  window size, in k-mers
  L  span of a window = k + w - 1   (the app's own "k-mer": kmc/kraken2/sshash k)

Density columns are named `expected_density_factor_{mean,std}` and, per
dataset, `<dataset>_particular_density_factor_{mean,std}`.
"""

ORDERS_ID_FIELDS = ("order", "k", "w", "L", "gm_k", "projection", "gm_table", "gm_lut", "mod_t")
ORDERS_DESC_FIELDS = ("description", "cli")
EXPECTED_FIELDS = ("expected_density_factor_mean", "expected_density_factor_std")


def span_length(k, w):
    """L = k + w - 1."""
    return int(k) + int(w) - 1


def particular_fields(datasets):
    out = []
    for ds in datasets:
        out += [f"{ds}_particular_density_factor_mean", f"{ds}_particular_density_factor_std"]
    return tuple(out)


def orders_fields(datasets, with_gm_lut=False, extra_density=()):
    """Column order of the orders sheet."""
    ids = tuple(f for f in ORDERS_ID_FIELDS if with_gm_lut or f != "gm_lut")
    return ids + ORDERS_DESC_FIELDS + EXPECTED_FIELDS + particular_fields(datasets) + tuple(extra_density)


def with_klw(row, k, w):
    """Set k, w and L = k + w - 1 on `row` (returns it)."""
    row["k"], row["w"], row["L"] = int(k), int(w), span_length(k, w)
    return row


def missing_cells(rows, fields, optional=()):
    """[(row_index, field)] of empty cells; `optional` fields (e.g. std of a
    deterministic order) are not reported. Used by each exporter to refuse
    writing an orders sheet with holes."""
    skip = set(optional)
    return [(i, f) for i, r in enumerate(rows) for f in fields
            if f not in skip and r.get(f) in (None, "")]
