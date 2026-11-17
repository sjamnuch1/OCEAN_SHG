#!/usr/bin/env python3
"""Compare link-space and real-space raw complex(DP) R X caches."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("link_cache", type=Path)
    parser.add_argument("realspace_cache", type=Path)
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    args = parser.parse_args()

    link = np.fromfile(args.link_cache, dtype=np.complex128)
    realspace = np.fromfile(args.realspace_cache, dtype=np.complex128)
    if link.size != realspace.size:
        parser.error(
            f"cache lengths differ: {link.size} complex values versus "
            f"{realspace.size}"
        )
    if link.size == 0:
        parser.error("the caches are empty")

    difference = realspace - link
    norm_link = np.linalg.norm(link)
    norm_realspace = np.linalg.norm(realspace)
    norm_difference = np.linalg.norm(difference)
    scale = max(norm_link, norm_realspace, np.finfo(float).tiny)
    relative = norm_difference / scale
    maximum = np.max(np.abs(difference))

    print(f"complex_values       {link.size}")
    print(f"link_norm            {norm_link:.16e}")
    print(f"realspace_norm       {norm_realspace:.16e}")
    print(f"difference_norm      {norm_difference:.16e}")
    print(f"relative_difference  {relative:.16e}")
    print(f"maximum_abs_error    {maximum:.16e}")
    print(f"tolerance            {args.rtol:.16e}")
    print("status               " + ("PASS" if relative <= args.rtol else "FAIL"))
    return 0 if relative <= args.rtol else 1


if __name__ == "__main__":
    raise SystemExit(main())

