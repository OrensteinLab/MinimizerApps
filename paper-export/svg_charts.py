"""SVG bar charts of the paper's results figures (make_figures.py): one chart per metric and dataset, the schemes as bars
from 0, the standard deviation over the seeds as a red error bar, the value above each bar, the change from the default
scheme inside the bar (or above a short one) and the significance stars above the labels.

Standard library only; the charts are plain SVG elements, composed into a figure by make_figures.py.
"""

import math

SANS = "DejaVu Sans, Helvetica, Arial, sans-serif"
# grey for the pseudo-random orders (default, xor, xor-pcg-hash), colours for the advanced schemes: GreedyMini green, the
# two mod-minimizers a light and a dark purple
COL = {"default": "#8a9099", "xor": "#d0d4d9", "xor-pcg-hash": "#a7adb4", "GreedyMini": "#1a9850",
       "mod-minimizer (xor)": "#9e8ad6", "mod-minimizer (GreedyMini)": "#5b3a9a"}
TICK = 22           # font size of the y-tick numbers
AXW = 74            # width left of a chart for the y-tick numbers
TICK_GAP = 32       # least distance between two y ticks (px), so that their numbers do not overlap
HEAD = 60           # default room above the highest tick for the value and change labels and the significance stars
ERR = "#c0392b"     # the standard deviation over the seeds, as a red error bar
STAGE2_OPACITY = 0.5   # stage 2 of a stacked bar: the colour of the scheme, lighter


class Fig:
    """A list of SVG elements of width W."""

    def __init__(self, W):
        self.W = W
        self.el = []

    def t(self, x, y, txt, size=18, anchor="middle", w="normal", fill="#222", rot=None, halo=False):
        """Text; halo: a white outlined copy underneath, so that a grid line behind the text does not cross it."""
        txt = txt.replace("&", "&amp;")
        tr = f' transform="rotate({rot} {x} {y})"' if rot else ""
        if halo:
            self.el.append(f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}" font-weight="{w}" fill="#fff" '
                           f'stroke="#fff" stroke-width="4" stroke-linejoin="round" font-family="{SANS}"{tr}>{txt}</text>')
        self.el.append(f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}" font-weight="{w}" fill="{fill}" '
                       f'font-family="{SANS}"{tr}>{txt}</text>')

    def r(self, x, y, w, h, fill, op=None):
        """A rectangle; op: fill opacity."""
        o = f' fill-opacity="{op}"' if op else ""
        self.el.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(w, 0):.1f}" height="{max(h, 0):.1f}" fill="{fill}" '
                       f'stroke="none" stroke-width="1" rx="0"{o}/>')

    def l(self, x1, y1, x2, y2, c="#999", sw=1.2):
        """A line."""
        self.el.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{c}" stroke-width="{sw}"/>')


def xlabel(f, x, base, lines, size=21):
    """A scheme name (its lines) under the bar at x, rotated by 45 degrees; every line ends below the bar, the next line
    XLABEL_STEP lower, i.e. one line height further across the text direction."""
    for j, ln in enumerate(lines):
        f.t(round(x + 4, 1), round(base + 12 + j * XLABEL_STEP * size, 1), ln, size, "end", "normal", "#444", rot=-45)


XLABEL_STEP = 1.1 * 2 ** 0.5   # the vertical step between the lines of a scheme name, per unit of font size


def fmtv(v, extra=0):
    """A value with 0 to 3 decimals by its magnitude, plus extra."""
    d = 0 if abs(v) >= 100 else (1 if abs(v) >= 10 else (2 if abs(v) >= 1 else 3))
    return f"{v:.{d + extra}f}"


def chart_fmt(values, ref, base=fmtv):
    """The value labels of one chart: base, with up to 3 more decimals for the whole chart where a value whose change from
    the default (ref) shows as nonzero would otherwise print the same as the default."""
    for e in range(4):
        if all(base(v, e) != base(ref, e) for v in values if ref and round(abs(v / ref - 1) * 100, 2) > 0):
            break
    return lambda v: base(v, e)


def ticks_for(lo, hi, zero, nmax=8):
    """(bottom, top, step) of a y axis for values in [lo, hi]: the finest step of 1, 2 or 5 x 10^e with at most nmax
    intervals; top is the first tick at or above hi. The axis starts at 0 (zero=True) or at the last tick below lo."""
    if zero:
        lo = 0.0
    mag = 10 ** math.floor(math.log10(max(hi - lo, 1e-12)))
    for st in sorted(k * mag * e for k in (1, 2, 5) for e in (0.01, 0.1, 1, 10)):
        st = float(f"{st:.12g}")
        top = float(f"{math.ceil(hi / st - 1e-9) * st:.12g}")
        below = math.floor(lo / st + 1e-9) * st
        bot = 0.0 if zero else float(f"{max(0.0, below - (st if below >= lo - 1e-12 else 0)):.12g}")
        if (top - bot) / st <= (nmax if zero else 7) + 1e-9:
            return bot, top, st


def adaptive(vals, px):
    """The axis from 0 to the first tick at or above the largest value (with its std), with the finest step of ticks_for
    that leaves TICK_GAP between the ticks of an axis px high."""
    return ticks_for(0.0, max(v + sd for _, v, sd, *_ in vals), True, max(2, int(px // TICK_GAP)))


def axis(f, x0, x1, base, sc, bot, top, st):
    """The y axis: grid lines and tick numbers from bot to top in steps of st."""
    f.l(x0, base, x1, base, "#999", 1)
    f.t(x0 - 5, base + 5, f"{bot:g}", TICK, "end", "normal", "#555")
    t = round(bot + st, 12)
    while t <= top * (1 + 1e-9) + 1e-15:
        yy = base - sc(t)
        f.l(x0, yy, x1, yy, "#e3e6ea", 1)
        f.t(x0 - 5, yy + 5, f"{t:.12g}", TICK, "end", "normal", "#555")
        t = round(t + st, 12)


def errbar(f, xc, ytop, sh):
    """A red error bar of half-height sh around the bar top ytop."""
    if not sh:
        return
    f.l(xc, ytop - sh, xc, ytop + sh, ERR, 2.2)
    f.l(xc - 6, ytop - sh, xc + 6, ytop - sh, ERR, 2.2)
    f.l(xc - 6, ytop + sh, xc + 6, ytop + sh, ERR, 2.2)


def lum(c, op=1.0):
    """Relative luminance of a #rrggbb fill drawn with opacity op on white."""
    r, g, b = (int(c[i:i + 2], 16) / 255 * op + (1 - op) for i in (1, 3, 5))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def pct_label(f, xc, y, v, ref, bw, inside):
    """The change from the default scheme's value ref, in % (two decimals), after an arrow pointing down (lower than the
    default) or up (higher). inside = (top of the bar, its height, its fill, opacity, half the error bar): drawn inside the
    bar, below the error bar, white or dark by the fill; above the bar (at y) when the bar is too short.
    Returns True when drawn above the bar."""
    if not ref:
        return False
    d = (v / ref - 1) * 100
    size = min(18, round(bw / 5.2, 1))
    txt = f"{abs(d):.2f}%"
    col, halo, above = "#444", True, True
    if inside[1] >= inside[4] + size * 1.9:
        y = inside[0] + inside[4] + size * 0.75 + 6
        col = "#fff" if lum(inside[2], inside[3]) < 0.55 else "#222"
        halo, above = False, False
    tw = len(txt) * size * 0.6                       # text width (estimate)
    aw = size * 0.75                                 # arrow width
    x0 = xc - (tw + aw + 3) / 2                      # the arrow and the text, centred together on xc
    xt, ax = x0 + aw + 3 + tw / 2, x0 + aw / 2
    f.t(xt, y, txt, size, "middle", "normal", col, halo=halo)
    h = size * 0.9                                   # the arrow, centred on the text's x-height
    yt = y - size * 0.35 - h / 2
    yb = yt + h
    hw, hh = aw / 2, size * 0.45
    tip, base = (yb, yt) if d < 0 else (yt, yb)
    sg = 1 if d < 0 else -1
    f.el.append(f'<line x1="{ax}" y1="{base}" x2="{ax}" y2="{tip - sg * hh * 0.6}" stroke="{col}" stroke-width="{size * 0.2:.1f}"/>')
    f.el.append(f'<path d="M{ax - hw},{tip - sg * hh} L{ax + hw},{tip - sg * hh} L{ax},{tip} z" fill="{col}"/>')
    return above


def stars_label(f, xc, y, stars):
    """The significance stars of a bar's change from the default scheme, above its labels (y: baseline of the topmost one)."""
    if stars:
        f.t(xc, y - 19, stars, 22, "middle", "bold", "#222", halo=True)


def minibars(f, x, y, w, h, title, vals, head=HEAD, xlabels=None):
    """Bars of one chart: vals = list of (scheme, value[, std[, stars]]); head: the room above the top tick; xlabels:
    {scheme: lines of its name} drawn under the axis, or None."""
    vals = [(t[0], t[1], (t[2] if len(t) > 2 else None) or 0.0, t[3] if len(t) > 3 else "") for t in vals if t[1] is not None]
    if title:
        f.t(x + AXW + (w - AXW) / 2, y + 20, title, 22, "middle", "bold")   # centred on the plot area
    base = y + h
    top = y + (32 if title else 14)
    bot, ymax, st = adaptive(vals, base - top - head)
    x0 = x + AXW
    k = (base - top - head) / (ymax - bot)
    sc = lambda v: k * (v - bot)
    ref = next((v for s, v, _, _ in vals if s == "default"), None)
    fmt = chart_fmt([u for _, u, _, _ in vals], ref)
    axis(f, x0, x + w, base, sc, bot, ymax, st)
    bw = (w - AXW) / (len(vals) + 0.6)
    gx = x0 + bw * 0.3
    for i, (s, v, sd, stars) in enumerate(vals):
        hh, eh, xc = sc(v), k * sd, gx + i * bw + bw / 2
        f.r(gx + i * bw + 3, base - hh, bw - 6, hh, COL[s])
        errbar(f, xc, base - hh, eh)
        f.t(xc, base - hh - eh - 5, fmt(v), 20, "middle", "normal", "#222", halo=True)
        if s != "default":
            above = pct_label(f, xc, base - hh - eh - 28, v, ref, bw, (base - hh, hh, COL[s], 1.0, eh))
            stars_label(f, xc, base - hh - eh - (28 if above else 5), stars)
        if xlabels:
            xlabel(f, xc, base, xlabels[s])
    f.l(x0, base, x + w, base)


def stacked_mini(f, x, y, w, h, title, vals, head=HEAD, xlabels=None):
    """Stacked bars of one chart (stage 1 dark, stage 2 light): vals = list of (scheme, stage 1, stage 2[, std of the
    total[, stars]]); the axis starts at 0, so that both stages stay visible."""
    vals = [(t[0], t[1], t[2], (t[3] if len(t) > 3 else None) or 0.0, t[4] if len(t) > 4 else "") for t in vals]
    if title:
        f.t(x + AXW + (w - AXW) / 2, y + 20, title, 22, "middle", "bold")
    base = y + h
    top = y + (32 if title else 14)
    x0 = x + AXW
    bot, ymax, st = ticks_for(0.0, max(a + b + sd for _, a, b, sd, _ in vals), True, max(2, int((base - top - head) // TICK_GAP)))
    sc = lambda v: (base - top - head) * v / ymax
    ref = next((a + b for s, a, b, _, _ in vals if s == "default"), None)
    fmt = chart_fmt([u + t for _, u, t, _, _ in vals], ref, lambda v, e: f"{v:.{e}f}")
    axis(f, x0, x + w, base, sc, bot, ymax, st)
    bw = (w - AXW) / (len(vals) + 0.6)
    gx = x0 + bw * 0.3
    for i, (s, a, b, sd, stars) in enumerate(vals):
        h1, h2, eh, xc = sc(a), sc(b), sc(sd), gx + i * bw + bw / 2
        f.r(gx + i * bw + 3, base - h1, bw - 6, h1, COL[s])
        f.r(gx + i * bw + 3, base - h1 - h2, bw - 6, h2, COL[s], op=STAGE2_OPACITY)
        errbar(f, xc, base - h1 - h2, eh)
        f.t(xc, base - h1 - h2 - eh - 5, fmt(a + b), 20, "middle", halo=True)
        if s != "default":
            above = pct_label(f, xc, base - h1 - h2 - eh - 28, a + b, ref, bw, (base - h1 - h2, h2, COL[s], STAGE2_OPACITY, eh))
            stars_label(f, xc, base - h1 - h2 - eh - (28 if above else 5), stars)
        if xlabels:
            xlabel(f, xc, base, xlabels[s])
    f.l(x0, base, x + w, base)
