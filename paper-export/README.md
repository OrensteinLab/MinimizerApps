# Paper export

The scripts that compute the paper's results tables and figures from the result sets of the
benchmarks (`benchmark/<method>/README.md`). Every number, test and chart of the paper's results
comes from these scripts; the result sets provide only the measured data.

| Script | Writes |
|---|---|
| `make_tables.py` | the supplementary tables (`supplementary_tables/S<nn>_*.csv` and their list in `supplementary_tables/README.md`); with `--labels FILE`, their LaTeX labels; with `--density-table FILE`, the rows of Table 2 |
| `make_figures.py` | the results figures, `results_runtime` and `results_space` (PDF, and their SVG sources) |
| `common.py` | shared by both: the result sets, the schemes and datasets the paper reports, and their names |
| `svg_charts.py` | the bar charts of the figures, drawn as SVG |
| `method_figures/` | the figures of the Methods, which illustrate the schemes on one window: `center_rule.py` (the center rule), `greedymini_order.py` (the canonical GreedyMini order, on the real table `w5_k6_canonical.gm`) and `canonical_mod_minimizer.py` (the canonical mod-minimizer) |

Both take the result set of each method:

```
R="kraken2=KRAKEN2_DIR minimap2=MINIMAP2_DIR kmc=KMC_DIR sshash=SSHASH_DIR"   # the result set of each method
python3 paper-export/make_tables.py --results $R [--labels supplementary_labels.tex] [--density-table table2.tex]
python3 paper-export/make_figures.py --results $R [--out DIR] [--vector-out DIR]
```

The figures of the Methods do not depend on the result sets:

```
for f in center_rule greedymini_order canonical_mod_minimizer; do
  python3 paper-export/method_figures/$f.py [--out DIR] [--vector-out DIR]
done
```

They write the PDF to `--out` (default: the current folder) and, with `--vector-out`, keep its SVG source.
`greedymini_order.py` reads `orders/GreedyMini/w5_k6_canonical.gm`, which is built from its `.txt` rank list
(`orders/README.md`, "Regenerating a `.gm` table").

The paper's result sets are the full runs of the four methods together with their profiling
runs (`benchmark/<method>/experiments/profile_run.sh`). The scripts read the CSVs that the export
of each method writes to `<results>/csvs/` and, for Minimap2's mapping profile, the
`map_profile.log` of each run. Python 3.7 or later (with `benchmark/common`), standard library
only; `make_figures.py` and the scripts of `method_figures/` also need `rsvg-convert` (librsvg) for the PDF.

## What the paper reports

- The schemes: the default scheme of each method, `xor-pcg-hash` (and `xor` in the tables),
  GreedyMini, and the mod-minimizer with a `xor` or a GreedyMini order on t-mers. Every
  GreedyMini-based scheme is reported with the polarity of Phi that maps A,T to 0 (`at0-cg1`),
  chosen arbitrarily; the table `polarity_comparison` compares it with the other polarity.
- Kraken2 at minimum hit groups 2 (its default); Minimap2 with the occurrence threshold fixed to
  that of the default scheme (and, in `minimap2_build_mapping`, also derived from each index); SSHash with the
  GreedyMini order built for 15-mers; KMC3 with the binary GreedyMini table.

## Statistics

Computed in the scripts with `benchmark/common`:

- Significance tests (`significance.py`; the tables `significance_runtime` and `significance_space` and the stars of the figures): every run of a
  scheme is paired with the run of the default scheme of the same dataset, setting and seed. Per
  scheme, dataset and metric: the geometric mean of the five ratios to the default scheme and its
  95% confidence interval, and a two-sided paired t-test on the logarithms of the ratios,
  Bonferroni-corrected over the schemes compared with the default scheme on the same dataset and
  metric. A change is significant if the corrected p-value is below 0.05; otherwise the
  performance is equivalent to that of the default scheme if the 90% confidence interval lies
  within +-5% (two one-sided tests), and inconclusive if not.
- The two polarities of Phi (`polarity_comparison`): the same test, pairing the run of each GreedyMini-based
  scheme with A,T mapped to 1 with its run with A,T mapped to 0 of the same seed, Bonferroni-
  corrected over the GreedyMini-based schemes of a dataset and metric.
- Correlations with density (`density_effect.py`; the tables `*_density_correlations`): the Pearson
  correlation of the ratios of the density factor and of a metric to those of the default
  scheme, over the paired runs of all datasets and of each dataset.
