#!/usr/bin/env python3
"""Check structural compatibility of the CV and core conduction spaces."""

from __future__ import annotations

import argparse
import os
from pathlib import Path


def parse_manifest(path: Path) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        words = raw.split()
        if words:
            result[words[0]] = words[1:]
    required = {
        "status",
        "sector",
        "completion",
        "dimensions",
        "input_brange",
        "kmesh",
    }
    missing = sorted(required - result.keys())
    if missing:
        raise ValueError(f"{path}: missing keys: {', '.join(missing)}")
    if result["status"] != ["RUNNING"] or result["sector"] != ["CORE_FROM_VALENCE"]:
        raise ValueError(f"{path}: not a CORE_FROM_VALENCE manifest")
    if result["completion"] != ["RBXC_CV_COMPLETE"]:
        raise ValueError(f"{path}: calculation is incomplete")
    return result


def read_ints(path: Path, count: int) -> list[int]:
    words = path.read_text(encoding="utf-8").split()
    if len(words) < count:
        raise ValueError(f"{path}: expected at least {count} integers")
    try:
        return [int(word) for word in words[:count]]
    except ValueError as exc:
        raise ValueError(f"{path}: invalid integer input") from exc


def regular_file(path: Path) -> None:
    if not path.is_file():
        raise ValueError(f"not a regular file: {path}")
    if path.stat().st_size <= 0:
        raise ValueError(f"empty file: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check valence/core conduction-space structural compatibility"
    )
    parser.add_argument("--cv-manifest", required=True, type=Path)
    parser.add_argument("--valence-con", required=True, type=Path)
    parser.add_argument("--core-con", required=True, type=Path)
    parser.add_argument("--core-brange", required=True, type=Path)
    parser.add_argument("--core-kmesh", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    if args.out.exists():
        raise ValueError(f"refusing to overwrite: {args.out}")
    metadata = parse_manifest(args.cv_manifest)
    dimensions = metadata["dimensions"]
    if dimensions[:3] != ["ncon", "nk", "nalpha"] or len(dimensions) != 6:
        raise ValueError("malformed dimensions record")
    ncon, nk, _ = map(int, dimensions[3:])
    cv_brange = list(map(int, metadata["input_brange"]))
    cv_kmesh = list(map(int, metadata["kmesh"]))
    if len(cv_brange) != 4 or len(cv_kmesh) != 3:
        raise ValueError("malformed input_brange or kmesh record")

    core_brange = read_ints(args.core_brange, 4)
    core_kmesh = read_ints(args.core_kmesh, 3)
    if cv_brange[2:] != core_brange[2:]:
        raise ValueError(
            f"conduction band ranges differ: CV={cv_brange[2:]} core={core_brange[2:]}"
        )
    if ncon != core_brange[3] - core_brange[2] + 1:
        raise ValueError("ncon is inconsistent with the core conduction band range")
    if cv_kmesh != core_kmesh or nk != core_kmesh[0] * core_kmesh[1] * core_kmesh[2]:
        raise ValueError("CV and core k meshes differ or nk is inconsistent")

    regular_file(args.valence_con)
    regular_file(args.core_con)
    if args.valence_con.stat().st_size != args.core_con.stat().st_size:
        raise ValueError("con.u2.dat byte sizes differ")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    fd = os.open(args.out, flags, 0o644)
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        handle.write("status BASIS_MATCH_COMPLETE\n")
        handle.write("method STRUCTURAL_METADATA_MATCH\n")
        handle.write("content_hash NOT_CHECKED\n")
        handle.write("bloch_gauge NOT_VERIFIED\n")
        handle.write(f"valence_bytes {args.valence_con.stat().st_size}\n")
        handle.write(f"core_bytes {args.core_con.stat().st_size}\n")
        handle.write(f"ncon {ncon}\n")
        handle.write(f"nk {nk}\n")
        handle.write("kmesh " + " ".join(map(str, cv_kmesh)) + "\n")
        handle.write("conduction_brange " + " ".join(map(str, cv_brange[2:])) + "\n")
        handle.write(f"cv_manifest {args.cv_manifest.resolve()}\n")
        handle.write(f"valence_con_file {args.valence_con.resolve()}\n")
        handle.write(f"core_con_file {args.core_con.resolve()}\n")
    print(f"Wrote structural conduction-space certificate: {args.out}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        raise SystemExit(f"FATAL: {exc}") from exc
