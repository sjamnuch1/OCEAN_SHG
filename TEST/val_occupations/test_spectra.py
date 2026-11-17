#!/usr/bin/env python3
"""Exercise the compiled ocean.x against the analytic two-pole response and a rescaled exchange kernel."""

from math import cos, exp, isclose, pi, sin, sqrt
from pathlib import Path
import struct
import subprocess
from tempfile import TemporaryDirectory
import sys

INPUTS = {
    'xmesh.ipt': '1 1 1',
    'kmesh.ipt': '1 1 2',
    'nspin': '1',
    'ZNL': '0 0 0',
    'nbuse.ipt': '2',
    'qinunitsofbvectors.ipt': '0.001 0 0',
    'brange.ipt': '1 1\n2 3',
    'avecsinbohr.ipt': '1 0 0\n0 1 0\n0 0 1',
    'cks.normal': '.false.',
    'mode': '1',
    'nelectron': '2',
    'epsilon': '1',
    'nval.h': '1',
    'tmel_selector': '1',
    'enk_selector': '0',
    'bloch_selector': '3',
    'cnbse.write_rhs': '.false.',
    'runlist': '1\n0 0 0 __ __ 0 0 VAL',
    'lflag': '0',
    'bflag': '0',
    'bande': '1',
    'semitda': '0',
    'bse.in': '0\n0\nhay\n2 101 0 18 0.1 0',
    'tmels.info': '1 2 3 2',
    'enkfile': '-0.02\n0.02 0.06\n-0.02\n0.02 0.06',
    'val_occupation.ipt': '3000 fixed 0',
    'val_occupation_energies.dat': '4 2 1\n0.5 -0.01 0.01 0.03 1.0\n0.5 -0.01 0.01 0.03 1.0',
}

HARTREE_EV = 27.211386245988
KT = 3000 * 8.617333262145e-5
EDGES = [-0.01 * HARTREE_EV, 0.01 * HARTREE_EV, 0.03 * HARTREE_EV]


def run(command, settings=None, bad_spectrum=False, bubble=False, scales=(1.0, 1.0), nelec=2):
    with TemporaryDirectory() as tmp:
        directory = Path(tmp)
        inputs = INPUTS.copy()
        inputs["nelectron"] = str(nelec)
        if settings is None:
            del inputs["val_occupation.ipt"]
            del inputs["val_occupation_energies.dat"]
        else:
            inputs["val_occupation.ipt"] = settings
        if bad_spectrum:
            inputs["val_occupation_energies.dat"] = inputs["val_occupation_energies.dat"].replace("0.03", "0.04")
        if bubble:
            inputs.update({"bflag": "1", "xmesh.ipt": "3 3 3", "avecsinbohr.ipt": "10 0 0\n0 10 0\n0 0 10"})
            # Plane waves with a nonzero G component produce a nonzero
            # exchange for each pair; identical shapes mix the two transitions.
            val = [component for _ in range(54) for component in (1.0, 0.0)]
            con = [
                component
                for _ in range(2)  # k points
                for scale in scales
                for x in [0, 1, 2] * 9
                for component in (scale * cos(2*pi*x/3), scale * sin(2*pi*x/3))
            ]
            (directory / "val.u2.dat").write_bytes(struct.pack(f"={len(val)}d", *val))
            (directory / "con.u2.dat").write_bytes(struct.pack(f"={len(con)}d", *con))
        for filename, contents in inputs.items():
            (directory / filename).write_text(contents + "\n")
        (directory / "ptmels.dat").write_bytes(struct.pack("=8d", *(scales[0], 0.0, scales[1], 0.0) * 2))
        result = subprocess.run(command, cwd=directory, capture_output=True, text=True, timeout=40)
        if bad_spectrum:
            assert result.returncode != 0, result.stdout
            assert "counting spectrum disagrees" in result.stdout, result.stdout + result.stderr
            return
        if result.returncode != 0 or not (directory / "opcons").exists():
            raise RuntimeError(result.stdout + result.stderr)
        rows = (directory / "opcons").read_text().splitlines()[1:]
        return [(float(row.split()[0]), float(row.split()[2])) for row in rows]


def f(energy, mu):
    return 1 / (1 + exp((energy - mu) / KT))


def check_ratio(legacy, thermal, mu):
    edge = EDGES
    assert len(legacy) == len(thermal) == 101
    for (omega, unweighted), (omega_thermal, weighted) in zip(legacy, thermal):
        assert omega == omega_thermal
        poles = [1 / ((omega - (edge[i] - edge[0])) ** 2 + 0.1 ** 2) for i in (1, 2)]
        expected = sum(poles[i] * (f(edge[0], mu) - f(edge[i + 1], mu)) for i in range(2)) / sum(poles)
        assert isclose(weighted / unweighted, expected, rel_tol=1e-6), (omega, weighted / unweighted, expected)


def main():
    executable = Path(sys.argv[1]).resolve()
    command = [*sys.argv[2:], str(executable)]
    legacy = run(command)
    check_ratio(legacy, run(command, "3000 fixed 0"), 0)
    # The full four-band spectrum sets mu even though the fourth band is
    # outside the optical transition space. The target intentionally is fractional.
    energies = [e * HARTREE_EV for e in (-0.01, 0.01, 0.03, 1.0)]
    lo, hi = -1.0, 1.0
    for _ in range(60):
        mu = (lo + hi) / 2
        if 2 * sum(f(e, mu) for e in energies) > 2.125:
            hi = mu
        else:
            lo = mu
    check_ratio(legacy, run(command, "3000 conserve 0", nelec=2.125), mu)
    run(command, "3000 conserve 0", bad_spectrum=True)
    print("Legacy, fixed-mu and conserving IPA spectra passed. Mismatched energies rejected.")
    # Scaling conduction coefficients independently scales each exchange-kernel
    # row and column. Also scale the bare source by hand for this reference run.
    scales = tuple(sqrt(f(EDGES[0], 0) - f(ec, 0)) for ec in EDGES[1:])
    unweighted = run(command, bubble=True)
    expected = run(command, bubble=True, scales=scales)
    weighted = run(command, "3000 fixed 0", bubble=True)
    assert any(abs(a-b) > 0.01 * abs(b) for (_, a), (_, b) in zip(expected, unweighted))
    for (omega, a), (_, b) in zip(weighted, expected):
        assert isclose(a, b, rel_tol=1e-7), (omega, a, b)
    print("Occupation-weighted exchange spectrum agrees with explicit row and column scaling.")


if __name__ == "__main__":
    main()
