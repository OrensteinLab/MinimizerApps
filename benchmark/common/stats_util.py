"""Scipy-free t-test helpers shared by the benchmarks' statistics (minimap2 analysis/stats.py,
common/order_tests.py): Student-t tail via the regularized incomplete beta function (Numerical
Recipes betacf), Welch's t-test and the Holm step-down correction."""
from __future__ import annotations

import math
import statistics


# -- t distribution (regularized incomplete beta, Numerical Recipes betacf) --

def _betacf(a: float, b: float, x: float) -> float:
    qab, qap, qam = a + b, a + 1, a - 1
    c, d = 1.0, 1 - qab * x / qap
    d = 1 / (d if abs(d) > 1e-300 else 1e-300)
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d; d = 1 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1 + aa / c if abs(c) > 1e-300 else 1e300
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d; d = 1 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1 + aa / c if abs(c) > 1e-300 else 1e300
        de = d * c
        h *= de
        if abs(de - 1) < 1e-12:
            break
    return h


def _betai(a: float, b: float, x: float) -> float:
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x))
    if x < (a + 1) / (a + b + 2):
        return bt * _betacf(a, b, x) / a
    return 1 - bt * _betacf(b, a, 1 - x) / b


def t_two_sided_p(t: float, df: float) -> float:
    return _betai(df / 2, 0.5, df / (df + t * t))


def t_quantile(p: float, df: float) -> float:
    """Upper quantile: t such that P(T > t) = 1 - p, by bisection."""
    lo, hi = 0.0, 1000.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if 1 - t_two_sided_p(mid, df) / 2 < p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def welch(x: list[float], y: list[float]) -> tuple[float, float, float, float]:
    """(mean diff x-y, se, df, two-sided p)."""
    mx, my = statistics.mean(x), statistics.mean(y)
    vx, vy = statistics.variance(x), statistics.variance(y)
    se2 = vx / len(x) + vy / len(y)
    if se2 == 0:
        return mx - my, 0.0, float("inf"), 0.0 if mx != my else 1.0
    df = se2 ** 2 / ((vx / len(x)) ** 2 / (len(x) - 1) + (vy / len(y)) ** 2 / (len(y) - 1))
    t = (mx - my) / math.sqrt(se2)
    return mx - my, math.sqrt(se2), df, t_two_sided_p(abs(t), df)


def holm(ps: list[float]) -> list[float]:
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    adj = [0.0] * len(ps)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(ps) - rank) * ps[i]))
        adj[i] = running
    return adj
