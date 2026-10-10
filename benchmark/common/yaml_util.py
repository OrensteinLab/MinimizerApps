"""Minimal YAML loader for experiment files (stdlib-only fallback)."""

from __future__ import annotations

import re
from typing import Any


def load_yaml(text: str) -> Any:
    try:
        import yaml  # type: ignore
    except ImportError:
        return _parse_yaml(text)
    return yaml.safe_load(text)


def _parse_yaml(text: str) -> Any:
    lines: list[tuple[int, str]] = []
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indent = len(raw) - len(raw.lstrip(" "))
        lines.append((indent, raw.strip()))

    if not lines:
        return {}

    def parse_block(i: int, indent: int) -> tuple[Any, int]:
        if i >= len(lines):
            return {}, i
        cur_indent, cur = lines[i]
        if cur.startswith("- "):
            return _parse_list(i, cur_indent)
        return _parse_mapping(i, cur_indent)

    def _parse_list(i: int, indent: int) -> tuple[list[Any], int]:
        out: list[Any] = []
        while i < len(lines):
            cur_indent, cur = lines[i]
            if cur_indent != indent or not cur.startswith("- "):
                break
            item_text = cur[2:].strip()
            if item_text:
                if ":" in item_text:
                    key, val = item_text.split(":", 1)
                    item = {key.strip(): _coerce(val.strip())}
                    i += 1
                    while i < len(lines) and lines[i][0] > indent:
                        sub, i = parse_block(i, lines[i][0])
                        if isinstance(item, dict) and isinstance(sub, dict):
                            item.update(sub)
                        else:
                            break
                    out.append(item)
                    continue
                out.append(_coerce(item_text))
                i += 1
                continue
            i += 1
            if i < len(lines) and lines[i][0] > indent:
                sub, i = parse_block(i, lines[i][0])
                out.append(sub)
            else:
                out.append(None)
        return out, i

    def _parse_mapping(i: int, indent: int) -> tuple[dict[str, Any], int]:
        out: dict[str, Any] = {}
        while i < len(lines):
            cur_indent, cur = lines[i]
            if cur_indent < indent:
                break
            if cur_indent > indent:
                break
            if cur.startswith("- "):
                break
            key, _, rest = cur.partition(":")
            key = key.strip()
            rest = rest.strip()
            if not rest:
                i += 1
                if i < len(lines) and lines[i][0] > cur_indent:
                    if lines[i][1].startswith("- "):
                        val, i = _parse_list(i, lines[i][0])
                    else:
                        val, i = _parse_mapping(i, lines[i][0])
                    out[key] = val
                else:
                    out[key] = {}
                continue
            out[key] = _coerce(rest)
            i += 1
        return out, i

    root, _ = parse_block(0, lines[0][0])
    return root if isinstance(root, dict) else {}


def _coerce(raw: str) -> Any:
    if not raw:
        return ""
    if raw in ("true", "True", "yes", "Yes"):
        return True
    if raw in ("false", "False", "no", "No"):
        return False
    if raw.startswith('"') and raw.endswith('"'):
        return raw[1:-1]
    if raw.startswith("'") and raw.endswith("'"):
        return raw[1:-1]
    if re.fullmatch(r"-?\d+", raw):
        return int(raw)
    if re.fullmatch(r"-?\d+\.\d+", raw):
        return float(raw)
    return raw
