"""Derived minimizer / super-k-mer metrics for benchmark export."""

from __future__ import annotations

from scripts.config import default_signature_len


def minimizer_window(k: int, signature_len: int | None = None) -> int:
    p = signature_len if signature_len is not None else default_signature_len()
    return int(k) - p + 1


def mean_read_length(
    k: int,
    n_total_kmers: int | float | None,
    n_reads: int | float | None,
) -> float | None:
    """K-mer-weighted mean read length: n_kmers/n_reads + (k-1)."""
    if n_total_kmers is None or n_reads is None:
        return None
    n_rd = int(n_reads)
    n_tot = int(n_total_kmers)
    if n_rd <= 0 or n_tot <= 0:
        return None
    return round(n_tot / n_rd + (int(k) - 1), 4)


def super_kmer_metrics(
    *,
    k: int,
    signature_len: int | None,
    n_total_kmers: int | float | None,
    n_super_kmers: int | float | None,
    n_reads: int | float | None = None,
) -> dict[str, float | int | None]:
    """Return density metrics or empty values when inputs are missing."""
    empty = {
        "minimizer_window_w": None,
        "avg_super_kmer_len": None,
        "avg_read_len": mean_read_length(k, n_total_kmers, n_reads),
        "minimizer_density": None,
        "density_factor": None,
    }
    if n_total_kmers is None or n_super_kmers is None:
        return empty
    n_tot = int(n_total_kmers)
    n_sup = int(n_super_kmers)
    if n_sup <= 0 or n_tot <= 0:
        empty["minimizer_window_w"] = minimizer_window(k, signature_len)
        return empty
    w = minimizer_window(k, signature_len)
    avg_len = n_tot / n_sup
    density = n_sup / n_tot
    return {
        "minimizer_window_w": w,
        "avg_super_kmer_len": round(avg_len, 4),
        "avg_read_len": mean_read_length(k, n_tot, n_reads),
        "minimizer_density": round(density, 6),
        "density_factor": round(density * (w + 1), 4),
    }


SUPER_KMER_FIELD_NAMES = (
    "minimizer_window_w",
    "avg_super_kmer_len",
    "avg_read_len",
    "minimizer_density",
    "density_factor",
)


def enrich_super_kmer_metrics(row: dict) -> dict:
    """Merge super-k-mer metrics into a run row (in place + return)."""
    k = row.get("k")
    if k is None:
        return row
    p = row.get("signature_len")
    if p is not None:
        p = int(p)
    metrics = super_kmer_metrics(
        k=int(k),
        signature_len=p,
        n_total_kmers=row.get("n_total_kmers"),
        n_super_kmers=row.get("n_super_kmers"),
        n_reads=row.get("n_reads"),
    )
    row.update(metrics)
    return row
