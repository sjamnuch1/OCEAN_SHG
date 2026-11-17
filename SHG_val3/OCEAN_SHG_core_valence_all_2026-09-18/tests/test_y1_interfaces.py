#!/usr/bin/env python3
"""Self-contained interface tests; add --contract-exe for the Fortran dot product."""

from __future__ import annotations

import argparse
import math
import struct
import subprocess
import sys
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parents[1]


def run_basis_test() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        manifest = root / "shg_rbxc_cv.info"
        manifest.write_text(
            "\n".join(
                [
                    "status RUNNING",
                    "sector CORE_FROM_VALENCE",
                    "operator R_CV=-conjg(xes)_hole",
                    "binary_layout raw_complex128_Fortran_ncon_nk_nalpha_no_padding",
                    "units bohr_squared_per_hartree",
                    "electron_frame unchanged_from_valence_k_plus_q",
                    "middle_axis z",
                    "right_axis z",
                    "dimensions ncon nk nalpha 2 4 4",
                    "input_brange 1 3 4 5",
                    "kmesh 2 2 1",
                    "ZNL 1 0 0",
                    "site C 1",
                    "eta_hartree_q_reduced 0.01 0 0 0",
                    "energy_count 1",
                    "RbXc 1 1.0 27.2114 2.0 rhs.bin",
                    "completion RBXC_CV_COMPLETE",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        payload = b"identical synthetic conduction basis\x00\x01"
        valence = root / "valence.con.u2.dat"
        core = root / "core.con.u2.dat"
        valence.write_bytes(payload)
        core.write_bytes(payload)
        brange = root / "brange.ipt"
        brange.write_text("1 3\n4 5\n", encoding="utf-8")
        kmesh = root / "kmesh.ipt"
        kmesh.write_text("2 2 1\n", encoding="utf-8")
        cert = root / "basis.info"
        command = [
            sys.executable,
            str(HERE / "verify_cv_core_basis.py"),
            "--cv-manifest",
            str(manifest),
            "--valence-con",
            str(valence),
            "--core-con",
            str(core),
            "--core-brange",
            str(brange),
            "--core-kmesh",
            str(kmesh),
            "--out",
            str(cert),
        ]
        subprocess.run(command, check=True, capture_output=True, text=True)
        text = cert.read_text(encoding="utf-8")
        assert "status BASIS_MATCH_COMPLETE" in text
        assert "conduction_brange 4 5" in text

        core.write_bytes(payload[:-1] + bytes([payload[-1] ^ 1]))
        mismatch = root / "mismatch.info"
        failed = subprocess.run(
            command[:-1] + [str(mismatch)], capture_output=True, text=True
        )
        assert failed.returncode != 0
        assert "hashes differ" in failed.stderr


def run_perl_stager_dry_run() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        for name in ("Common", "DFT", "OPF", "PREP/BSE", "SCREEN"):
            (root / name).mkdir(parents=True, exist_ok=True)
        common = {
            "calc": {"mode": "xas", "edges": ["1 1 0"]},
            "bse": {"core": {"gmres": {"elist": [100.0]}}},
        }
        import json

        (root / "Common/postDefaultsOceanDatafile").write_text(
            json.dumps(common), encoding="utf-8"
        )
        payload = b"same-conduction-frame"
        (root / "PREP/BSE/con.u2.dat").write_bytes(payload)
        valence = root / "valence-con.u2.dat"
        valence.write_bytes(payload)
        vector = root / "cv.bin"
        vector.write_bytes(b"synthetic-vector")
        manifest = root / "cv.info"
        manifest.write_text(
            "status RUNNING\nsector CORE_FROM_VALENCE\n"
            "middle_axis z\nright_axis z\nenergy_count 1\n"
            "RbXc 1 1.0 27.2114 2.0 cv.bin\n"
            "completion RBXC_CV_COMPLETE\n",
            encoding="utf-8",
        )
        result = subprocess.run(
            [
                "perl",
                str(HERE / "prepare_y1_cv.pl"),
                "--ocean-run-root",
                str(root),
                "--target",
                str(root / "GV_FRESH"),
                "--component",
                "xzz",
                "--cv-manifest",
                str(manifest),
                "--valence-con",
                str(valence),
                "--dry-run",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert "previous" not in result.stdout.lower()
        assert "fresh GV stage is valid" in result.stdout


def write_complex_stream(path: Path, values: list[complex]) -> None:
    with path.open("wb") as handle:
        for value in values:
            handle.write(struct.pack("=dd", value.real, value.imag))


def run_contraction_test(executable: Path) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        d = [1 + 2j, -0.5 + 0.25j]
        y = [3 - 1j, 2 + 4j]
        expected = sum(a.conjugate() * b for a, b in zip(d, y))
        dfile, yfile = root / "d.bin", root / "y.bin"
        write_complex_stream(dfile, d)
        write_complex_stream(yfile, y)
        dnorm = math.sqrt(sum(abs(x) ** 2 for x in d))
        ynorm = math.sqrt(sum(abs(x) ** 2 for x in y))
        common = "dimensions ncon nk nalpha 2 1 1\nkmesh 1 1 1\ncore_edge_ZNL 1 0 0\n"
        dm = root / "d.info"
        dm.write_text(
            "status XC_MINUS_COMPLETE\nsector CORE_XAS\nbinary_layout "
            "raw_complex_DP_stream_in_OCEAN_psi_vtor_order\nright_axis x\n"
            "core_site C 1 K 1\n"
            + common
            + f"dipole {dnorm:.17e} {dfile}\n",
            encoding="utf-8",
        )
        ym = root / "y.info"
        ym.write_text(
            "status Y1_GV_COMPLETE\nsector CORE_FROM_VALENCE\nbinary_layout "
            "raw_complex_DP_stream_in_OCEAN_psi_vtor_order\ncomponent xyz\n"
            "left_axis x\nmiddle_axis y\nright_axis z\ncore_site C 1 K 1\n"
            + common
            + f"y1_norm {ynorm:.17e}\noutput_vector {yfile}\n"
            "source_index 1\nsource_omega_hartree 1\nsource_omega_ev 27.2114\n"
            "target_energy_hartree 2\ntarget_energy_ev 54.4228\n"
            "target_is_two_omega T\neta_hartree 0.01\n",
            encoding="utf-8",
        )
        output = root / "T1.dat"
        subprocess.run(
            [str(executable), str(dm), str(dfile), str(ym), str(yfile), str(output)],
            check=True,
            capture_output=True,
            text=True,
        )
        fields = {
            words[0]: words[1]
            for line in output.read_text(encoding="utf-8").splitlines()
            if len(words := line.split()) == 2
        }
        got = complex(float(fields["T1_real"]), float(fields["T1_imag"]))
        assert abs(got - expected) < 1.0e-13


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract-exe", type=Path)
    args = parser.parse_args()
    run_basis_test()
    run_perl_stager_dry_run()
    if args.contract_exe:
        run_contraction_test(args.contract_exe.resolve())
    print("All available interface tests passed")


if __name__ == "__main__":
    main()
