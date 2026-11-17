#!/usr/bin/env python3
"""Algebraic link/FFT equivalence and common photon-q phase cancellation."""

from __future__ import annotations

import numpy as np


def linear_index(index: np.ndarray, mesh: np.ndarray) -> int:
    return int((index[0] * mesh[1] + index[1]) * mesh[2] + index[2])


def evaluate_case(phase_shift: np.ndarray):
    rng = np.random.default_rng(4271)
    kmesh = np.array([3, 4, 5])
    xmesh = np.array([2, 3, 2])
    nk = int(np.prod(kmesh))
    nx = int(np.prod(xmesh))
    ncon = 3
    nalpha = 2
    tau = np.array([0.17, 0.31, 0.08])
    k0 = np.array([0.5, 0.0, 0.5])
    avec = np.array([[7.1, 1.2, 0.4], [0.3, 6.4, 0.8], [0.2, 0.5, 8.0]])
    cart_axis = 1
    weight = kmesh * avec[cart_axis, :] / (4.0 * np.pi)

    # New con.u2.dat order: primitive x fastest, then y, then z.
    fraction = np.array(
        [
            np.array([ix, iy, iz]) / xmesh
            for iz in range(xmesh[2])
            for iy in range(xmesh[1])
            for ix in range(xmesh[0])
        ]
    )
    # OCEAN vector order: k3 fastest, then k2, then k1.
    kindex = np.array(
        [
            [i1, i2, i3]
            for i1 in range(kmesh[0])
            for i2 in range(kmesh[1])
            for i3 in range(kmesh[2])
        ]
    )

    state = np.empty((nx, ncon, nk), dtype=np.complex128)
    for ik in range(nk):
        qmat, _ = np.linalg.qr(
            rng.normal(size=(nx, ncon)) + 1j * rng.normal(size=(nx, ncon))
        )
        state[:, :, ik] = qmat
    source = rng.normal(size=(ncon, nk, nalpha)) + 1j * rng.normal(
        size=(ncon, nk, nalpha)
    )

    full_bloch = np.empty_like(state)
    for ik, mesh_index in enumerate(kindex):
        phase = np.exp(
            1j
            * 2.0
            * np.pi
            * np.sum(((mesh_index + k0) / kmesh + phase_shift) * fraction, axis=1)
        )
        full_bloch[:, :, ik] = state[:, :, ik] * phase[:, None]

    link_result = np.zeros_like(source)
    for ialpha in range(nalpha):
        for ik, mesh_index in enumerate(kindex):
            for direction in range(3):
                for step in (-1, 1):
                    neighbor = mesh_index.copy()
                    neighbor[direction] = (
                        neighbor[direction] + step
                    ) % kmesh[direction]
                    ik_neighbor = linear_index(neighbor, kmesh)
                    transport_phase = np.exp(
                        -1j
                        * step
                        * 2.0
                        * np.pi
                        * fraction[:, direction]
                        / kmesh[direction]
                    )
                    link = np.conj(full_bloch[:, :, ik]).T @ (
                        transport_phase[:, None]
                        * full_bloch[:, :, ik_neighbor]
                    )
                    core_phase = np.exp(
                        1j
                        * step
                        * 2.0
                        * np.pi
                        * tau[direction]
                        / kmesh[direction]
                    )
                    link_result[:, ik, ialpha] += (
                        1j
                        * step
                        * weight[direction]
                        * core_phase
                        * (link @ source[:, ik_neighbor, ialpha])
                    )

    realspace_result = np.zeros_like(source)
    for ialpha in range(nalpha):
        for ix, xi in enumerate(fraction):
            phase = np.exp(
                1j
                * 2.0
                * np.pi
                * np.sum(((kindex + k0) / kmesh + phase_shift) * xi, axis=1)
            )
            contracted = np.einsum(
                "ck,ck->k", state[ix, :, :], source[:, :, ialpha]
            )
            field_r = np.fft.ifftn((phase * contracted).reshape(*kmesh))

            rho_fd = np.empty(tuple(kmesh), dtype=float)
            for n1 in range(kmesh[0]):
                for n2 in range(kmesh[1]):
                    for n3 in range(kmesh[2]):
                        cell = np.array([n1, n2, n3])
                        rho_fd[n1, n2, n3] = np.sum(
                            2.0
                            * weight
                            * np.sin(
                                2.0
                                * np.pi
                                * (cell + xi - tau)
                                / kmesh
                            )
                        )

            field_k = np.fft.fftn(rho_fd * field_r).reshape(-1)
            realspace_result[:, :, ialpha] += (
                np.conj(state[ix, :, :] * phase[None, :]) * field_k[None, :]
            )

    difference = realspace_result - link_result
    relative = np.linalg.norm(difference) / np.linalg.norm(link_result)
    maximum = np.max(np.abs(difference))
    expectation = np.vdot(source, realspace_result)
    relative_imaginary_expectation = abs(expectation.imag) / max(
        abs(expectation), np.finfo(float).tiny
    )

    return link_result, realspace_result, relative, maximum, relative_imaginary_expectation


def main() -> int:
    # These are uniform phase offsets, NOT substitutions of DFT states from
    # another mesh. A common exp(i*q.r) factor cancels for this local operator.
    # The actual core code uses the base-mesh phase; the nonzero cases also
    # check that removing the old common-q factor preserves its algebra.
    phase_shifts = [np.zeros(3), np.array([0.013, -0.017, 0.021]), np.array([0.2, 0.1, -0.3])]
    baseline = None
    passed = True
    for shift in phase_shifts:
        link, realspace, relative, maximum, imaginary = evaluate_case(shift)
        if baseline is None:
            baseline = (link.copy(), realspace.copy())
        link_change = np.linalg.norm(link-baseline[0]) / np.linalg.norm(baseline[0])
        realspace_change = np.linalg.norm(realspace-baseline[1]) / np.linalg.norm(baseline[1])
        print(f"phase_shift                    {shift}")
        print(f"relative_difference             {relative:.16e}")
        print(f"maximum_abs_error              {maximum:.16e}")
        print(f"relative_imaginary_expectation  {imaginary:.16e}")
        print(f"link_change_from_zero_phase     {link_change:.16e}")
        print(f"realspace_change_from_zero_phase {realspace_change:.16e}")
        passed = passed and max(relative, imaginary, link_change, realspace_change) < 1.0e-13
    print("status                          " + ("PASS" if passed else "FAIL"))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())

