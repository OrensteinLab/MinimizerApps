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


def write_runs(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        w.writerows(rows)


def test_points_pair_by_setting_and_seed(tmp_path):
    csvs = tmp_path / "csvs"
    csvs.mkdir()
    rows = []
    # the default's time differs by seed and by dataset, so a pairing across either gives other ratios
    for ds, scale in (("cod", 1), ("kestrel", 10)):
        for seed, t in (("0", 100), ("1", 200)):
            rows.append({"dataset": ds, "order": "native", "seed_idx": seed, "density_factor": "2.0",
                         "build_wall_s": str(t * scale)})
            rows.append({"dataset": ds, "order": "greedymini", "gm_k": "15", "seed_idx": seed,
                         "density_factor": "1.5", "build_wall_s": str(0.8 * t * scale)})
    write_runs(csvs / "sshash_runs.csv", rows)
    p = dict(de.points(tmp_path, "sshash"))["build_wall_s"]
    assert len(p) == 4
    assert all(x["metric_ratio"] == 0.8 and x["density_ratio"] == 0.75 for x in p)
    assert sorted((x["setting"], x["seed_idx"]) for x in p) == [("cod", "0"), ("cod", "1"), ("kestrel", "0"), ("kestrel", "1")]


def test_minimap2_pairs_only_the_fixed_threshold_block(tmp_path):
    csvs = tmp_path / "csvs"
    csvs.mkdir()
    run = lambda order, thr, t, d: {"dataset": "chrX", "k": "19", "w": "10", "mid_occ_threshold": thr,
                                    "order": order, "seed_idx": "0", "density_factor": d, "map_run_s": t}
    rows = [run("native", "159", "10", "2.0"),          # default-filter block
            run("xor", "159", "50", "2.0"),              # derived threshold that happens to equal native's
            run("native", "159", "12", "2.0"),           # fixed block: native's second run of the seed
            run("greedymini", "159", "9.6", "1.6")]
    write_runs(csvs / "minimap2_runs.csv", rows)
    p = dict(de.points(tmp_path, "minimap2"))["map_run_s"]
    assert [x["scheme"] for x in p] == ["greedymini"]
    assert p[0]["metric_ratio"] == 0.8 and p[0]["density_ratio"] == 0.8
