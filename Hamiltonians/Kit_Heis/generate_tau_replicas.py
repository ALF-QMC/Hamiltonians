#!/usr/bin/env python3
"""Generate replicas for bootstrap error estimates from covariance matrix eigen-decomposition.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def parse_input_file(file_path: Path) -> tuple[int, int, float, int, str, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Read input file with format used by generate_tau_replicas.F90."""
    with file_path.open("r", encoding="utf-8") as f:
        header_tokens = f.readline().split()
        if len(header_tokens) < 5:
            raise ValueError("Header must contain: Lt nbins beta Norb Channel")

        lt = int(header_tokens[0])
        nbins = int(header_tokens[1])
        beta = float(header_tokens[2])
        norb = int(header_tokens[3])
        channel = header_tokens[4]

        tau = np.empty(lt, dtype=np.float64)
        gt = np.empty(lt, dtype=np.float64)
        gt_err = np.empty(lt, dtype=np.float64)

        for i in range(lt):
            row = f.readline().split()
            if len(row) < 3:
                raise ValueError(f"Missing tau/gt/gt_err data at row {i + 1}")
            tau[i] = float(row[0])
            gt[i] = float(row[1])
            gt_err[i] = float(row[2])

        cov_tokens: list[str] = []
        for line in f:
            parts = line.split()
            if parts:
                cov_tokens.extend(parts)

    required = lt * lt
    if len(cov_tokens) < required:
        raise ValueError(f"Covariance matrix has {len(cov_tokens)} values, expected at least {required}")

    cov_values = np.array(cov_tokens[:required], dtype=np.float64)
    covtt = cov_values.reshape((lt, lt))

    return lt, nbins, beta, norb, channel, tau, gt, gt_err, covtt


def write_replica(
    out_path: Path,
    lt: int,
    nbins: int,
    beta: float,
    norb: int,
    channel: str,
    tau: np.ndarray,
    gt_new: np.ndarray,
    gt_err: np.ndarray,
) -> None:
    """Write one replica file with the same format as the Fortran program."""
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", encoding="utf-8") as f:
        f.write(f"{lt} {nbins} {beta} {norb} {channel}\n")
        for i in range(lt):
            f.write(f"{tau[i]:14.7f}  {gt_new[i]:16.10f}  {gt_err[i]:16.10f}\n")




def generate_tau_replicas(input_file: Path, number_of_replicas: int) -> None:
    """Generate replicas for bootstrap error estimates from covariance matrix eigen-decomposition."""

    if number_of_replicas < 1:
        raise ValueError("number_of_replicas must be >= 1")

    print(f"Reading file: {input_file}")
    lt, nbins, beta, norb, channel, tau, gt, gt_err, covtt = parse_input_file(input_file)
    print(lt, nbins, beta, norb, channel)

    en, u = np.linalg.eigh(covtt)
    sig = np.sqrt(en)

    # rang_wrap() in Fortran draws from N(0, 1).
    rng = np.random.default_rng()

    input_name = input_file.name
    for n_eta in range(1, number_of_replicas + 1):
        eta = rng.normal(loc=0.0, scale=sig, size=lt)
        gt_new = gt + u @ eta

        replica_dir = input_file.parent / f"{input_name}_{n_eta:03d}"
        file_out = replica_dir / "g_dat"
        print(f"Writing file: {file_out}")
        write_replica(file_out, lt, nbins, beta, norb, channel, tau, gt_new, gt_err)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate replicas for bootstrap error estimates "
            "from covariance matrix eigen-decomposition."
    )
    parser.add_argument("input_file", type=Path, help="Input file (e.g. Spintot)")
    parser.add_argument("number_of_replicas", type=int, help="Number of replicas to generate")
    args = parser.parse_args()
    generate_tau_replicas(args.input_file, args.number_of_replicas)

if __name__ == "__main__":
    main()
