"""Generic YAML multi-config experiment-sweep engine, shared by
benchmark/kraken2/experiments.py and benchmark/kmc/scripts/experiments.py.

A YAML experiment file holds a `defaults:` mapping plus a `configs:` list of
overrides; each config entry, merged over `defaults`, becomes one sweep run.
This module only knows how to load/merge/list those files and render a
config as a `VAR=val; ...` shell-env line -- each app supplies its own
`run_keys` (the config field names it understands) and `var_mapping` (how
those fields map to its own shell variable names), since the two apps'
benchmark harnesses take different parameters.
"""

from __future__ import annotations

import shlex
import sys
from pathlib import Path
from typing import Any, Mapping, NoReturn, Sequence

from yaml_util import load_yaml


def die(msg: str) -> NoReturn:
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(1)


def resolve_path(experiments_dir: Path, name_or_path: str) -> Path:
    p = Path(name_or_path)
    if p.is_file():
        return p.resolve()
    for candidate in (
        experiments_dir / name_or_path,
        experiments_dir / f"{name_or_path}.yaml",
        experiments_dir / f"{name_or_path}.yml",
    ):
        if candidate.is_file():
            return candidate.resolve()
    die(f"experiment not found: {name_or_path} (look in {experiments_dir})")


def load_experiment(path: Path) -> dict[str, Any]:
    data = load_yaml(path.read_text())
    if not isinstance(data, dict):
        die(f"{path}: expected mapping at top level")
    configs = data.get("configs")
    if not isinstance(configs, list) or not configs:
        die(f"{path}: missing non-empty configs list")
    return data


def _as_str(v: Any) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, (list, tuple)):
        return ",".join(str(x) for x in v)
    return str(v)


def expand_configs(spec: dict[str, Any], run_keys: Sequence[str]) -> list[dict[str, str]]:
    defaults = spec.get("defaults") or {}
    if not isinstance(defaults, dict):
        die("defaults must be a mapping")
    out: list[dict[str, str]] = []
    for i, entry in enumerate(spec["configs"]):
        if not isinstance(entry, dict):
            die(f"configs[{i}] must be a mapping")
        merged = {**defaults, **entry}
        out.append({k: _as_str(merged.get(k, "")) for k in run_keys})
    return out


def list_experiments(experiments_dir: Path) -> list[tuple[str, str, int]]:
    rows: list[tuple[str, str, int]] = []
    for path in sorted(experiments_dir.glob("*.yaml")):
        spec = load_experiment(path)
        desc = str(spec.get("description", "")).strip().replace("\n", " ")
        rows.append((path.stem, desc, len(spec["configs"])))
    return rows


def shell_env_line(
    cfg: dict[str, str],
    var_mapping: Mapping[str, str],
    bool_vars: Sequence[str] = (),
) -> str:
    parts: list[str] = []
    for src, dst in var_mapping.items():
        val = cfg.get(src, "")
        if dst in bool_vars:
            parts.append(f'{dst}={1 if val in ("1", "true", "True") else 0}')
        else:
            parts.append(f"{dst}={shlex.quote(val)}")
    return "; ".join(parts)


def shell_env_lines(
    spec: dict[str, Any],
    run_keys: Sequence[str],
    var_mapping: Mapping[str, str],
    bool_vars: Sequence[str] = (),
) -> list[str]:
    """One shell_env_line() per config, in one pass over the whole spec --
    lets a caller expand + render every config with a single subprocess
    instead of re-parsing the YAML once per config index."""
    return [
        shell_env_line(cfg, var_mapping, bool_vars)
        for cfg in expand_configs(spec, run_keys)
    ]
