"""Correlation footer rows shared by the summary sheets (xlsx_book `footer_rows`).

For every `{metric}_ratio` column of a summary sheet, one row under the table per
driver (e.g. the density ratio, the distinct-minimizer ratio) holding the Pearson
correlation of that column with the driver's column over all rows. The cell is
Excel's own CORREL (live and traceable) with the computed value cached, so a
viewer that does not recalculate still shows it. Ratios are this order / the
baseline order of the same setting, so rows of different settings are comparable.
"""
from __future__ import annotations

from typing import Sequence

from xlsx_book import DATA_START_ROW
from xlsx_util import col_letter


def pearson(xs: Sequence[float], ys: Sequence[float]) -> float | None:
    """Pearson r; None when undefined (fewer than 3 points or a constant series)."""
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 1e-18 or syy <= 1e-18:
        return None
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return round(max(-1.0, min(1.0, sxy / (sxx * syy) ** 0.5)), 3)


def _num(v):
    try:
        return None if v is None or v == "" or isinstance(v, dict) else float(v)
    except (TypeError, ValueError):
        return None


def correlation_rows(
    fields: Sequence[str],
    rows: Sequence[dict],
    drivers: Sequence[tuple[str, str]],
    columns: frozenset | None = None,
    value_of=None,
) -> list[dict]:
    """xlsx_book footer rows {"label", "cells": [(field, formula, cached)]}, one per
    (label, driver column) in `drivers` that is present in `fields`. `columns`
    restricts the `_ratio` columns that get a cell; `value_of(row, field)` supplies
    the number of a cell when rows hold formulas instead of values."""
    n = len(rows)
    if not n:
        return []
    get = value_of or (lambda r, f: r.get(f))
    first, last = DATA_START_ROW, DATA_START_ROW + n - 1
    fields = list(fields)
    out = []
    for label, driver in drivers:
        if driver not in fields:
            continue
        dl = col_letter(fields.index(driver) + 1)
        cells: list[tuple[str, str, float | None]] = []
        for c, f in enumerate(fields):
            if not f.endswith("_ratio") or (columns is not None and f not in columns):
                continue
            cl = col_letter(c + 1)
            pairs = [(_num(get(r, driver)), _num(get(r, f))) for r in rows]
            pairs = [(a, b) for a, b in pairs if a is not None and b is not None]
            value = pearson([a for a, _ in pairs], [b for _, b in pairs])
            # a constant column has no correlation: CORREL gives #DIV/0!, shown blank
            formula = f'IFERROR(CORREL(${dl}${first}:${dl}${last},{cl}${first}:{cl}${last}),"")'
            cells.append((f, formula, value))
        out.append({"label": label, "cells": cells})
    return out
