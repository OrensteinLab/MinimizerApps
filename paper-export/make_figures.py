#!/usr/bin/env python3
"""The results figures of the paper, computed from the result sets of the benchmarks.

    make_figures.py --results kraken2=DIR minimap2=DIR kmc=DIR sshash=DIR [--out DIR] [--vector-out DIR]

Writes --out/results_runtime.pdf and --out/results_space.pdf (the figures of the paper; default: paper-export/figures) and
their SVG sources to --vector-out (default: --out). Each figure has a panel per method (A Kraken2, B Minimap2,
C KMC3, D SSHash) and a bar chart per metric and dataset: the mean over the seeds of each scheme, its standard deviation
as a red error bar, and its change from the default scheme, marked *, ** or *** when significant (Bonferroni-adjusted
p < 0.05, 0.01, 0.001; benchmark/common/significance.py, computed here from the runs). Every GreedyMini-based scheme is
shown with the polarity at0-cg1 of the projection (A,T -> 0), chosen arbitrarily.

Each DIR is a result set of benchmark/<app> (the paper's: the full runs with the profiling runs,
benchmark/<app>/experiments/profile_run.sh); the figures read DIR/csvs, the CSVs that the export of each app writes.
Standard library only, plus rsvg-convert (librsvg) for the PDF.
"""

import argparse
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path[:0] = [HERE, os.path.join(REPO, "benchmark", "common")]
import significance  # noqa: E402
import svg_charts as C  # noqa: E402
from common import APPS, ResultSets, fl, key_label, label, minimap2_blocks  # noqa: E402

SCH4 = ["default", "xor-pcg-hash", "GreedyMini", "mod-minimizer (xor)", "mod-minimizer (GreedyMini)"]
SCM = ["default", "xor-pcg-hash", "GreedyMini", "mod-minimizer (xor)"]   # Minimap2 at k = 15 (ONT): no GreedyMini O_t
SCK = ["default", "xor-pcg-hash", "GreedyMini"]                          # KMC3: no mod-minimizer
XLAB = {"default": ["default"], "xor-pcg-hash": ["xor-pcg-hash"], "GreedyMini": ["C-GreedyMini"],
        "mod-minimizer (xor)": ["C-mod (xor)"], "mod-minimizer (GreedyMini)": ["C-mod (GM)"]}   # short, so that each
# name stays under its bar
# Minimap2: (dataset, chart title key, t of the mod-minimizers, schemes); the mod-minimizer with a GreedyMini order on
# t-mers only at k = 19 (t = 9, the w20_k9 order)
MM2_SETS = (("chm13_chrX", "chm13_chrx", 9, SCH4), ("grch38_wg", "grch38", 9, SCH4), ("chm13_rel3_ont", "ont", 5, SCM))
KMC_DS = ["f_vesca", "g_gallus", "m_balbisiana", "h_sapiens_2", "h_sapiens_3"]
# SSHash: three of the six collections, each showing a different effect on the query (few lookups; many lookups; heavy
# buckets); all six are in the supplementary tables. (dataset, chart title key)
SSH_DS = [("ncbi-virus", "ncbi viruses"), ("kestrel", "kestrel"), ("hprc", "hprc")]

H = 230             # height of a chart's plot area (px)
YL = 58             # width left of the y-tick numbers for the y-axis label (up to two lines)
TITLE_GAP = 0       # room (px) left between the highest label of a row and the titles above it
TITLE = {"library": "strain-exclusion library", "cami": "CAMI II", "fda_argos": "FDA-ARGOS",
         "chm13_chrx": "CLR, CHM13 chrX (<i>k</i> = 19)", "grch38": "CLR, GRCh38 (<i>k</i> = 19)", "ont": "ONT, GRCh38 (<i>k</i> = 15)",
         "f_vesca": "<i>F. vesca</i>", "g_gallus": "<i>G. gallus</i>", "m_balbisiana": "<i>M. balbisiana</i>",
         "h_sapiens_2": "<i>H. sapiens</i> 2", "h_sapiens_3": "<i>H. sapiens</i> 3",
         "shared": " ",     # the room of a title, drawn once over the row (space)
         "kestrel": "kestrel", "hprc": "HPRC", "ncbi viruses": "NCBI viruses"}
YLAB = {"build_time": "Build time (s)", "classify_processing": "Classification|time (s)", "map_processing": "Mapping|time (s)",
        "query_processing": "Streaming|query time (s)", "count_time": "Runtime (s)", "database_size": "Database size (GB)",
        "build_memory": "Build peak|memory (GB)", "index_size": "Index size (GB)", "map_memory": "Mapping peak|memory (GB)",
        "tmp_disk": "Temporary disk|space (GB)", "bits_per_kmer": "Bits per L-mer"}
GREY = 'font-size="22" font-weight="normal" fill="#555"'
KMC_LEGEND = (f'<rect x="0" y="-18" width="20" height="20" fill="#777"/><text x="27" y="0" font-size="22" fill="#555" '
              f'font-family="{C.SANS}">stage 1 (distribution)</text><rect x="290" y="-18" width="20" height="20" fill="#777" '
              f'fill-opacity="{C.STAGE2_OPACITY}"/><text x="317" y="0" font-size="22" fill="#555" font-family="{C.SANS}">stage 2 (sorting and counting)</text>')


# ---------------------------------------------------------------- data

def by_scheme(rows, app, keep=lambda x: True):
    """{scheme: {metric: (mean, ratio to the default, std)}} of the summary rows kept, from their <metric>_mean, _ratio and
    _std columns; the first row of a scheme."""
    d = {}
    for x in rows:
        if not keep(x):
            continue
        s = label(x, app)
        if s is None or s in d:
            continue
        d[s] = {k[:-5]: (fl(x[k]), fl(x.get(k[:-5] + "_ratio", "")), fl(x.get(k[:-5] + "_std", "")))
                for k in x if k.endswith("_mean") and fl(x[k]) is not None}
    return d


class Data:
    """The summaries of the four result sets and the significance of each scheme's change from the default."""

    def __init__(self, R):
        kr = R.csv("kraken2", "kraken2_summary.csv")
        self.kraken2_build = by_scheme(R.csv("kraken2", "kraken2_build.csv"), "kraken2")
        self.kraken2_cls = {ds: by_scheme(kr, "kraken2", lambda x, ds=ds: x["dataset"] == ds and x["minimum_hit_groups"] == "2")
                            for ds in ("cami", "fda-argos")}
        mm = minimap2_blocks(R.csv("minimap2", "minimap2_summary.csv"))[0]   # at the default scheme's occurrence threshold
        self.minimap2 = {ds: by_scheme(mm, "minimap2", lambda x, ds=ds: x["dataset"] == ds) for ds, *_ in MM2_SETS}
        km = R.csv("kmc", "kmc3_summary.csv")
        self.kmc = {ds: by_scheme(km, "kmc", lambda x, ds=ds: x["dataset"] == ds) for ds in KMC_DS}
        sh = R.csv("sshash", "sshash_summary.csv")
        self.sshash = {ds: by_scheme(sh, "sshash", lambda x, ds=ds: x["dataset"] == ds) for ds, _ in SSH_DS}
        self.sig = {app: significance.rows_of(str(R[app]), app) for app in APPS}

    def stars(self, app, setting, metric, s):
        """*, **, *** when the change of the scheme s from the default is significant (Bonferroni-adjusted p < 0.05, 0.01,
        0.001) on the setting (the first fields of significance.py's setting, e.g. "cami/2") and metric."""
        rows = [r for r in self.sig[app] if r["metric"] == metric and key_label(r["scheme"], app) == s
                and (r["setting"] == setting or r["setting"].startswith(setting + "/"))]
        assert len(rows) <= 1, ("several tests of one bar", app, setting, metric, s)
        if rows and rows[0].get("verdict") in ("lower", "higher"):
            p = float(rows[0]["p_bonferroni"])
            return "***" if p < 0.001 else "**" if p < 0.01 else "*"
        return ""


def bars(data, d, key, schemes, sig):
    """The bars of a chart: (scheme, mean, std, stars) of each scheme that has the metric key; sig = (app, setting)."""
    return [(s, d[s][key][0], d[s][key][2] or 0, data.stars(*sig, key, s)) for s in schemes if s in d and key in d[s]]


# ---------------------------------------------------------------- layout

def svg_text(t):
    """<i>...</i> in a title or label as an italic tspan."""
    return t.replace("<i>", '<tspan font-style="italic">').replace("</i>", "</tspan>")


def xlabels(t, schemes):
    """The scheme names on the x axis, with the t of the mod-minimizers once, as a second line of the rightmost one."""
    lab = {s: list(v) for s, v in XLAB.items()}
    mods = [s for s in schemes if s.startswith("mod-minimizer")]
    if mods:
        lab[mods[-1]].append(f'<tspan font-style="italic">t</tspan> = {t}')
    return lab


def xlab_h(lines_of, size=21):
    """The room below the axis for the scheme names of a chart (lines_of: their lines): the lowest point of the rotated
    labels, as C.xlabel draws them."""
    low = 0
    for lines in lines_of:
        for j, ln in enumerate(lines):
            low = max(low, 12 + j * C.XLABEL_STEP * size + 0.7071 * len(re.sub("<[^>]+>", "", ln)) * size * 0.52)   # 0.52: mean glyph width / size
    return round(low + 12)


def plot_px(o):
    """The plot area (base to top tick) of a chart with the default room C.HEAD above it."""
    return H + (24 if o["title"] else 0) - (32 if o["title"] else 14) - C.HEAD


def head_need(chart, W, P):
    """The room above the top tick that the labels of one chart take (value, change and stars above the bars, as drawn by
    C.minibars and C.stacked_mini), for a plot area P px high from the base to the top tick."""
    if chart and chart[0] == "stacked":
        vals = [(t[0], t[1] + t[2], (t[3] if len(t) > 3 else None) or 0.0, t[4] if len(t) > 4 else "", t[2]) for t in chart[1]]
        bot = 0.0
        _, ymax, _ = C.ticks_for(0.0, max(v + sd for _, v, sd, _, _ in vals), True, max(2, int(P // C.TICK_GAP)))
    else:
        vals = [(t[0], t[1], (t[2] if len(t) > 2 else None) or 0.0, t[3] if len(t) > 3 else "", None) for t in chart if t[1] is not None]
        bot, ymax, _ = C.adaptive([v[:4] for v in vals], P)
    k = P / (ymax - bot)
    bw = (W - C.AXW) / (len(vals) + 0.6)
    size = min(18, round(bw / 5.2, 1))
    ref = next((v for s, v, *_ in vals if s == "default"), None)
    need = 12                                                    # the number of the top tick
    for s, v, sd, st, top in vals:
        hh, eh, ext = k * (v - bot), k * sd, 20                  # ext: the value label
        if s != "default" and ref:
            above = not ((k * top if top is not None else hh) >= eh + size * 1.9)   # the change label above the bar
            ext = 28 + size * 0.8 if above else ext
            if st:
                ext = (28 if above else 5) + 19 + 16             # the stars above the labels
        need = max(need, ext - (P - hh - eh))
    return need


def opts(o):
    return {"ylab": True, "title": True, "xlab": True, **(o[0] if o else {})}


def row_head(charts, W):
    """The room above the top tick for a row of charts: the largest need of its charts plus TITLE_GAP, so that the titles
    sit close to their plots while the plot areas keep their size and stay aligned."""
    return round(max(head_need(c, W, plot_px(opts(o))) for _, _, c, *o in charts) + TITLE_GAP)


def make_charts(charts, W, lab, head=None):
    """The charts of one row, each as (width, height, SVG elements); charts: list of (metric, title key, chart[, opts]),
    a chart being a list of bars or ("stacked", list of stacked bars), opts a dict of ylab (the metric left of the axis),
    title (the dataset above the plot) and xlab (the scheme names, lab, below the axis), all True by default."""
    head = row_head(charts, W) if head is None else head
    out = []
    for metric, ds, c, *o in charts:
        o = opts(o)
        yl = YL if o["ylab"] else 0
        h = plot_px(o) + (32 if o["title"] else 14) + head
        f = C.Fig(W + 20 + yl)
        title = svg_text(TITLE[ds]) if o["title"] else ""
        stacked = c and c[0] == "stacked"
        draw = C.stacked_mini if stacked else C.minibars
        draw(f, 10 + yl, 0, W, h, title, c[1] if stacked else c, head=head, xlabels=lab if o["xlab"] else None)
        if o["ylab"]:                                # centered on the y axis (from its top tick to the base), by the ticks
            lines = YLAB[metric].replace("L-mer", "<i>L</i>-mer").split("|")
            yc = ((32 if title else 14) + head + h) / 2
            for j, ln in enumerate(lines):
                f.t(round(10 + yl - 6 - (len(lines) - 1 - j) * 24, 1), round(yc, 1), svg_text(ln), 21, "middle", "normal", "#333", rot=-90)
        names = [b[0] for b in (c[1] if stacked else c)]
        out.append((f.W, h + (xlab_h([lab[s] for s in names]) if o["xlab"] else 14), f.el))
    return out


def separated(chart, top=8):
    """A chart with a thin vertical rule at its left edge, separating it from a chart of another metric on its left."""
    w, h, el = chart
    return (w, h, [f'<line x1="2" y1="{top}" x2="2" y2="{h - 8}" stroke="#bbb" stroke-width="1.5"/>'] + el)


def grid(metrics, sets, chart):
    """Charts of a grid of metrics (rows) by datasets (columns): chart(metric, dataset)."""
    return [(m, ds, chart(m, ds), dict(ylab=j == 0, title=i == 0, xlab=i == len(metrics) - 1))
            for i, m in enumerate(metrics) for j, ds in enumerate(sets)]


def kw(k=None, w=None):
    """(k = .., w = ..) after a panel title, as SVG tspans (spaces at tspan boundaries collapse, so dx)."""
    if k is None:
        return (f'<tspan dx="16" {GREY}>(</tspan><tspan {GREY} font-style="italic">w</tspan><tspan dx="5" {GREY}>=</tspan>'
                f'<tspan dx="5" {GREY}>{w})</tspan>')
    return (f'<tspan dx="16" {GREY}>(</tspan><tspan {GREY} font-style="italic">k</tspan><tspan dx="5" {GREY}>=</tspan>'
            f'<tspan dx="5" {GREY}>{k},</tspan><tspan dx="8" {GREY} font-style="italic">w</tspan><tspan dx="5" {GREY}>=</tspan>'
            f'<tspan dx="5" {GREY}>{w})</tspan>')


def compose(panels):
    """One figure of panels as SVG: panels = [(title, legend, rows)], rows a list of rows of charts (make_charts), each
    panel headed by its letter (A, B, ...), a legend (an SVG fragment, or None) centered below its charts."""
    el, y, Wmax = [], 0, 0
    for n, (title, legend, rows) in enumerate(panels):
        i = title.find("<") % (len(title) + 1)                   # the panel letter, bold, before the method
        title = f'{"ABCDEFGH"[n]}<tspan dx="16">{title[:i]}</tspan>{title[i:]}'
        el.append(f'<text x="12" y="{y + 30}" font-size="28" font-weight="bold" fill="#222" font-family="{C.SANS}">{title}</text>')
        y += 44
        for row in rows:
            x = 0
            for w_, h_, e in row:
                el.append(f'<g transform="translate({x},{y})">' + "".join(e) + "</g>")
                x += w_
            Wmax = max(Wmax, x)
            y += max(h_ for _, h_, _ in row)
        if legend:
            el.append(f'<g transform="translate({x / 2 - 310},{y + 14})">{legend}</g>')
            y += 30
        y += 12
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{Wmax}" height="{y}" viewBox="0 0 {Wmax} {y}">'
            f'<rect width="100%" height="100%" fill="#fff"/>' + "".join(el) + "</svg>")


def minimap2_rows(data, metrics, W):
    """Minimap2 on the simulated reads (CHM13 chrX, GRCh38; k = 19) and the real ONT reads (k = 15), the occurrence
    threshold fixed to the default scheme's: a row per metric = (figure metric, column)."""
    out = []
    for i, (m, key) in enumerate(metrics):
        cs = [(m, title, bars(data, data.minimap2[ds], key, sch, ("minimap2", ds)),
               {"ylab": j == 0, "title": i == 0, "xlab": i == len(metrics) - 1})
              for j, (ds, title, t, sch) in enumerate(MM2_SETS)]
        hd = row_head(cs, W)                                     # one room above the plots for the row
        for (ds, title, t, sch), c in zip(MM2_SETS, cs):
            out += make_charts([c], W, xlabels(t, sch), head=hd)
    return [out[:3], out[3:]]


def sshash_rows(data, metrics):
    """SSHash: a grid of metrics (rows; (figure metric, column)) by the three collections."""
    col = dict(metrics)
    ds_of = {title: ds for ds, title in SSH_DS}
    return grid([m for m, _ in metrics], [title for _, title in SSH_DS],
                lambda m, title: bars(data, data.sshash[ds_of[title]], col[m], SCH4, ("sshash", ds_of[title])))


def runtime(data):
    """Figure results_runtime: Kraken2 build and classification (minimum hit groups 2; processing time, without the
    database load), Minimap2 build and mapping, KMC3 counting (two stages), SSHash build and streaming queries."""
    build = data.kraken2_build
    k2 = make_charts([("build_time", "library", bars(data, build, "build_wall_s", SCH4, ("kraken2", "strain_exclusion")))] +
                     [("classify_processing", title, bars(data, data.kraken2_cls[ds], "classify_run_s", SCH4, ("kraken2", ds + "/2")),
                       {"ylab": j == 0}) for j, (title, ds) in enumerate((("cami", "cami"), ("fda_argos", "fda-argos")))],
                     537, xlabels(6, SCH4))
    k2[1] = separated(k2[1])                                     # build | classification: two metrics in one row
    mm = minimap2_rows(data, (("build_time", "build_wall_s"), ("map_processing", "map_run_s")), 556)
    # the bar is the total (wall) time, which is tested; stage 1 dark, the rest (stage 2) light
    km = make_charts(grid(["count_time"], KMC_DS, lambda m, ds: ("stacked", [
        (s, data.kmc[ds][s]["count_stage1_s"][0], data.kmc[ds][s]["count_wall_s"][0] - data.kmc[ds][s]["count_stage1_s"][0],
         data.kmc[ds][s]["count_wall_s"][2],
         data.stars("kmc", ds, "count_wall_s", s)) for s in SCK if s in data.kmc[ds]])), 326, XLAB)
    sh = make_charts(sshash_rows(data, [("build_time", "build_wall_s"), ("query_processing", "query_run_s")]), 552, xlabels(10, SCH4))
    return compose([("Kraken2" + kw(31, 5), None, [k2]), ("Minimap2" + kw(w=10), None, mm),
                    ("KMC3" + kw(9, 20), KMC_LEGEND, [km]), ("SSHash" + kw(21, 11), None, [sh[:3], sh[3:]])])


def space(data):
    """Figure results_space: Kraken2 database size and build memory, Minimap2 index size and mapping memory, KMC3
    temporary disk, SSHash bits per L-mer and build memory."""
    build = data.kraken2_build
    k2 = make_charts([("database_size", "shared", bars(data, build, "index_size_gb", SCH4, ("kraken2", "strain_exclusion"))),
                      ("build_memory", "shared", bars(data, build, "build_max_rss_gb", SCH4, ("kraken2", "strain_exclusion")))],
                     816, xlabels(6, SCH4))
    w0, h0, e0 = k2[0]                                           # one title for both charts (one dataset), centered on the row
    k2[0] = (w0, h0, e0 + [f'<text x="{w0}" y="20" font-size="22" text-anchor="middle" font-weight="bold" fill="#222" '
                           f'font-family="{C.SANS}">{TITLE["library"]}</text>'])
    k2[1] = separated(k2[1], top=34)                             # database size | build peak memory, below the title
    mm = minimap2_rows(data, (("index_size", "index_size_gb"), ("map_memory", "map_max_rss_gb")), 556)
    km = make_charts(grid(["tmp_disk"], KMC_DS, lambda m, ds: bars(data, data.kmc[ds], "tmp_disk_gb", SCK, ("kmc", ds))),
                     326, XLAB)
    sh = make_charts(sshash_rows(data, [("bits_per_kmer", "bits_per_kmer"), ("build_memory", "build_max_rss_gb")]), 552, xlabels(10, SCH4))
    return compose([("Kraken2" + kw(31, 5), None, [k2]), ("Minimap2" + kw(w=10), None, mm),
                    ("KMC3" + kw(9, 20), None, [km]), ("SSHash" + kw(21, 11), None, [sh[:3], sh[3:]])])


def write(name, svg, out, vector_out):
    """<vector_out>/<name>.svg and <out>/<name>.pdf."""
    src = os.path.join(vector_out, name + ".svg")
    with open(src, "w") as fh:
        fh.write(svg)
    pdf = os.path.join(out, name + ".pdf")
    subprocess.check_call(["rsvg-convert", "-f", "pdf", src, "-o", pdf])
    print(name, "->", pdf)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", nargs=4, required=True, metavar="APP=DIR", help="the result set of each of " + ", ".join(APPS))
    ap.add_argument("--out", default=os.path.join(HERE, "figures"), help="folder of the PDFs (default: %(default)s)")
    ap.add_argument("--vector-out", default=None, help="folder of the SVG sources (default: --out)")
    a = ap.parse_args(argv)
    try:
        R = ResultSets(a.results)
    except ValueError as e:
        ap.error(str(e))
    vector_out = a.vector_out or a.out
    for d in (a.out, vector_out):
        os.makedirs(d, exist_ok=True)
    data = Data(R)
    write("results_runtime", runtime(data), a.out, vector_out)
    write("results_space", space(data), a.out, vector_out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
