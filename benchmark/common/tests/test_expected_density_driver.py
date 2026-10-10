import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import expected_density_driver as d  # noqa: E402

R = d.Rounding(sample=4, avg=4, std=4, single_std=None)


def test_aggregate():
    assert R.aggregate([]) == {"n_seeds": 0, "density_factor_mean": None, "density_factor_std": None}
    assert R.aggregate([2.0]) == {"n_seeds": 1, "density_factor_mean": 2.0, "density_factor_std": None}
    a = d.Rounding(4, 4, 4, 0.0).aggregate([1.0, 3.0])
    assert (a["density_factor_mean"], a["density_factor_std"]) == (2.0, 1.0)
    assert d.Rounding(4, 4, 4, 0.0).aggregate([2.0])["density_factor_std"] == 0.0


def test_compute_rows_unavailable():
    def prepare(v):
        if v == "missing":
            raise d.Unavailable
    def sample(v, ctx, seed):
        if seed == 2:
            raise d.Unavailable
        return float(seed)
    rows = d.compute_rows(["a", "missing"], seeds_of=lambda v: [1, 2, 3], sample=sample, prepare=prepare,
                          row_id=lambda v: {"order": v}, rounding=R)
    assert rows[0] == {"order": "a", "n_seeds": 2, "density_factor_mean": 2.0, "density_factor_std": 1.0}
    assert rows[1]["n_seeds"] == 0 and rows[1]["density_factor_mean"] is None


def test_random_dna_is_deterministic_and_keys_length():
    s = d.random_dna(100)
    assert s == d.random_dna(100) and set(s) <= {0, 1, 2, 3}
    keys = d.canonical_keys(s, 5, lambda v, seed: v ^ seed, lambda x, l: x, 7)
    assert len(keys) == 96


def test_ensure_cache_exact_and_merge(tmp_path):
    path = tmp_path / "c.csv"
    fields = ("order", *d.STAT_FIELDS)
    calls = []

    def compute(todo):
        calls.append(list(todo))
        return [{"order": v, "n_seeds": 1, "density_factor_mean": 1.5, "density_factor_std": None} for v in todo]

    def run(wanted, **kw):
        return d.ensure_cache(
            path, fields, wanted, key_of_variant=lambda v: v, key_of_row=lambda r: r.get("order") or None,
            load=lambda: d.read_csv_rows(path, lambda r: r.get("order") or None), compute=compute, **kw)

    out = run(["a", "b"])
    assert set(out) == {"a", "b"} and out["a"]["density_factor_std"] == ""
    run(["a", "b"])
    assert len(calls) == 1  # hit
    run(["a", "b", "c"])  # exact mode: different variant set -> recompute everything
    assert calls[-1] == ["a", "b", "c"]
    run(["a", "d"], merge=True)  # merge: only the missing one
    assert calls[-1] == ["d"]
    assert set(run(["a"], merge=True)) == {"a", "b", "c", "d"} and len(calls) == 3
    run(["a"], force=True, stamp={"order": "z"})
    assert set(d.read_csv_rows(path, lambda r: r["order"])) == {"z"}
