import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import migrate_scheme_names as m  # noqa: E402


def make_set(root):
    run = root / "runs" / "cod" / "mod-minimizer-xor" / "t10" / "0"
    run.mkdir(parents=True)
    (run / "meta.json").write_text(json.dumps({
        "dataset": "cod", "order": "mod-minimizer-xor", "order_seed": 7,
        "command": ["sshash build -i x -O mod-minimizer-xor --mod-t 10 -o y", "sshash query -i y -q z"]}, indent=2) + "\n")
    (run / "build.log").write_text("[M::main] CMD: minimap2 -d i.mmi --minimizer-order greedymini --order-seed 1 r.fa\nrest\n")
    csvs = root / "csvs"
    csvs.mkdir()
    (csvs / "app_orders.csv").write_text("order,k,n_orders\nxor,15,1\n")
    (csvs / "app_density_expected_orders.csv").write_text("order,density_factor_mean\nxor,2.0\n")
    (csvs / "app_runs.csv").write_text("dataset,order,signature_order\ncod,xor,custom\n")
    return run


def test_migrates_meta_logs_and_csvs(tmp_path):
    run = make_set(tmp_path)
    m.migrate("sshash", tmp_path)
    meta = json.loads((run / "meta.json").read_text())
    assert list(meta)[:3] == ["dataset", "scheme", "order_seed"] and meta["scheme"] == "mod-minimizer-xor"
    assert meta["command"][0] == "sshash build -i x --scheme mod-minimizer-xor --mod-t 10 -o y"
    csvs = tmp_path / "csvs"
    assert sorted(p.name for p in csvs.iterdir()) == ["app_density_expected_schemes.csv", "app_runs.csv", "app_schemes.csv"]
    assert (csvs / "app_schemes.csv").read_text().splitlines()[0] == "scheme,k,n_schemes"
    assert (csvs / "app_runs.csv").read_text().splitlines()[0] == "dataset,scheme,signature_order"
    m.migrate("minimap2", tmp_path)
    assert "--scheme greedymini --order-seed 1" in (run / "build.log").read_text()


def test_second_run_changes_nothing(tmp_path):
    make_set(tmp_path)
    m.migrate("minimap2", tmp_path)
    assert m.migrate("minimap2", tmp_path) == []


def test_minimap2_interim_flag_is_migrated_too(tmp_path):
    run = tmp_path / "runs" / "chrX" / "xor" / "0"
    run.mkdir(parents=True)
    (run / "build.log").write_text("[M::main] CMD: minimap2 -d i.mmi --minimizer-scheme xor --order-seed 1 "
                                   "--greedymini-file /x/common/minimizer-schemes/t.gm r.fa\n")
    m.migrate("minimap2", tmp_path)
    assert (run / "build.log").read_text() == ("[M::main] CMD: minimap2 -d i.mmi --scheme xor --order-seed 1 "
                                              "--greedymini-file /x/common/minimizer-schemes/t.gm r.fa\n")
