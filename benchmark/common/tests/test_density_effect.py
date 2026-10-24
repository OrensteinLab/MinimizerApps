"""density_effect.py: paired ratios to the default of the same setting and seed, and the correlation."""
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import density_effect as de  # noqa: E402


def test_correlation_overall_and_per_setting():
    pts = [{"scheme": "s%d" % (i % 3), "density_ratio": d, "metric_ratio": 1 + 0.8 * (d - 1)}
           for i, d in enumerate([0.7, 0.75, 0.8, 0.85, 0.9, 0.95])]
    for p in pts:
        p["setting"] = "cod" if p["density_ratio"] < 0.82 else "kestrel"
    rows = de.rows_of("sshash", "build_wall_s", pts)
    assert [r["setting"] for r in rows] == ["all", "cod", "kestrel"]
    assert rows[0]["r"] == 1.0 and rows[0]["schemes"] == 3 and rows[0]["n"] == 6
    assert [r["n"] for r in rows[1:]] == [3, 3] and all(r["r"] == 1.0 for r in rows[1:])


def test_points_pair_by_setting_and_seed(tmp_path):
    csvs = tmp_path / "csvs"
    csvs.mkdir()
    rows = []
    for seed in ("0", "1"):
        rows.append({"dataset": "cod", "order": "native", "seed_idx": seed, "density_factor": "2.0",
                     "build_wall_s": "100", "build_foreign_busy_cores": "0"})
        rows.append({"dataset": "cod", "order": "greedymini", "gm_k": "15", "seed_idx": seed,
                     "density_factor": "1.5", "build_wall_s": "80", "build_foreign_busy_cores": "2"})
    with open(csvs / "sshash_runs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        w.writerows(rows)
    got = dict(de.points(tmp_path, "sshash"))
    p = got["build_wall_s"]
    assert len(p) == 2 and p[0]["metric_ratio"] == 0.8 and p[0]["density_ratio"] == 0.75


def test_default_with_several_runs_uses_their_mean(tmp_path):
    csvs = tmp_path / "csvs"
    csvs.mkdir()
    rows = [{"dataset": "chrX", "mid_occ_threshold": "159", "k": "19", "order": "native", "seed_idx": "0",
             "density_factor": "2.0", "map_run_s": t} for t in ("10", "12")]
    rows.append({"dataset": "chrX", "mid_occ_threshold": "159", "k": "19", "order": "greedymini", "seed_idx": "0",
                 "density_factor": "1.6", "map_run_s": "9.9"})
    with open(csvs / "minimap2_runs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    p = dict(de.points(tmp_path, "minimap2"))["map_run_s"]
    assert len(p) == 1 and p[0]["metric_ratio"] == 0.9 and p[0]["density_ratio"] == 0.8
