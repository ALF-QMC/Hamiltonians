#!/usr/bin/env python3
"""Calculate the magnetotropic susceptibilty k composed out of tree terms k1, k2, k3."""

import os
import subprocess
from multiprocessing import Pool
from functools import partial
import math
from pathlib import Path
import argparse

import numpy as np
import yaml
import f90nml

from py_alf.ana import read_scal

from generate_tau_replicas import generate_tau_replicas

def calc_k2(dir_name="."):
    """Calculate second term of magnetotropic susceptibility k2 from the output of Max_SAC."""

    # Read parameters from file
    params = f90nml.read(f'{dir_name}/parameters')

    # Read g_dat file (first line has mixed types: int int float int [str])
    with open(f'{dir_name}/g_dat', 'r', encoding="utf-8") as f:
        header = f.readline().split()
        beta = float(header[2])
        data = np.loadtxt(f)

    # tau = data[:, 0]
    # Stau = data[:, 1]
    # Stau_err = data[:, 2]

    # Read Best_fit file
    best_fit = np.loadtxt(f'{dir_name}/Best_fit', encoding="utf-8")
    omi = best_fit[:, 0]
    ai = best_fit[:, 1]

    # Calculate k2
    k2 = 0.0
    for i in range(params['VAR_Max_Stoch']['Ngamma']):
        if omi[i] == 0.0:
            print(f"Warning: omi[{i}] is zero, skipping"
                "this term in k2 calculation to avoid division by zero.")
            print(omi[i], ai[i])
            continue  # Avoid division by zero; contribution is zero in this case.
        k2 += 2.0 * data[0, 1] * ai[i] / omi[i] * np.tanh(beta * omi[i] / 2.0)

    result = k2 / float(params['VAR_Lattice']['L1'] * params['VAR_Lattice']['L2'] * 2)
    with open(f"{dir_name}/k2.dat", "w", encoding="utf-8") as f:
        f.write(f"{result:.8E}\n")
    return result


def _jk_mean_and_error_real(values: np.ndarray) -> tuple[float, float]:
    """Return real jackknife mean and error from jackknife sample values."""
    n_bins = values.shape[0]
    if n_bins < 2:
        raise ValueError("Need at least 2 bins for jackknife error estimation.")

    mean = np.mean(values)
    err = math.sqrt((n_bins - 1) * np.mean((values - mean) ** 2))
    return mean, err


def get_n_m(b_qmc: float, theta: float) -> tuple[np.ndarray, np.ndarray]:
    """Calculate n and m vectors.

    Parameters:
    b_qmc: Magnetic field strength
    theta: Magnetic field angle in units of pi (0.5 means 90 degrees)
    
    n = e X B
    m = e X (e X B)
    """

    # e1 = 1.0 / math.sqrt(6.0)
    # e2 = 1.0 / math.sqrt(6.0)
    # e3 = -2.0 / math.sqrt(6.0)

    e1 = -1.0 / math.sqrt(2.0)
    e2 = 1.0 / math.sqrt(2.0)
    e3 = 0.0

    # e1 = 1.0 / math.sqrt(3.0)
    # e2 = 1.0 / math.sqrt(3.0)
    # e3 = 1.0 / math.sqrt(3.0)

    b1 = b_qmc * np.sin(theta * np.pi) / np.sqrt(6.0) \
        + b_qmc * np.cos(theta * np.pi) / np.sqrt(3.0)
    b2 = b_qmc * np.sin(theta * np.pi) / np.sqrt(6.0) \
        + b_qmc * np.cos(theta * np.pi) / np.sqrt(3.0)
    b3 = -2.0 * b_qmc * np.sin(theta * np.pi) / np.sqrt(6.0) \
        + b_qmc * np.cos(theta * np.pi) / np.sqrt(3.0)

    n = np.zeros(3, dtype=float)
    m = np.zeros(3, dtype=float)

    n[0] = -(e2 * b3 - e3 * b2)
    n[1] = -(e3 * b1 - e1 * b3)
    n[2] = -(e1 * b2 - e2 * b1)

    m[0] = e2 * (e1 * b2 - e2 * b1) - e3 * (e3 * b1 - e1 * b3)
    m[1] = e3 * (e2 * b3 - e3 * b2) - e1 * (e1 * b2 - e2 * b1)
    m[2] = e1 * (e3 * b1 - e1 * b3) - e2 * (e2 * b3 - e3 * b2)

    return n, m


def calc_k1k3(directory: Path) -> None:
    """Run the k1/k3 analysis and write k1k3.dat."""
    params = f90nml.read(f'{directory}/parameters')

    kit_heis = params["VAR_Kit_Heis"]

    volume = float(params["VAR_Lattice"]["L1"] * params["VAR_Lattice"]["L2"] * 2)

    n, m = get_n_m(float(kit_heis["Hab"]), float(kit_heis["Htheta"]))

    g1 = 2.3
    # g2 = 2.3
    g3 = 1.3

    g_a = (g3 + 2.0 * g1) / 3.0
    g_c = (g3 - g1) / 3.0

    obs_x, sign, _ = read_scal(directory, "Esx_scal")
    obs_y, sign, _ = read_scal(directory, "Esy_scal")
    obs_z, sign, _ = read_scal(directory, "Esz_scal")
    n_jacks = obs_x.shape[0]

    s = np.zeros((3, n_jacks), dtype=float)

    s[0] = g_a * obs_x[:, 0].real + g_c * (obs_y[:, 0].imag + obs_z[:, 0].real)
    s[1] = g_a * obs_y[:, 0].imag + g_c * (obs_x[:, 0].real + obs_z[:, 0].real)
    s[2] = g_a * obs_z[:, 0].real + g_c * (obs_x[:, 0].real + obs_y[:, 0].imag)

    k1_bins = np.zeros((n_jacks,), dtype=float)
    for i in range(3):
        k1_bins += m[i] * s[i]
    k1, k1_err = _jk_mean_and_error_real(k1_bins / sign)

    k3_bins = np.zeros((n_jacks,), dtype=float)
    for i in range(3):
        for j in range(3):
            k3_bins += s[i] * s[j] * float(kit_heis["Beta"]) * n[i] * n[j]
    k3, k3_err = _jk_mean_and_error_real(k3_bins / sign**2)

    return k1 / volume, k1_err / volume, k3 / volume, k3_err / volume


def _run_max_sac(alf_dir, sim_dir, dir_name):
    print(f"Running Max_SAC in {dir_name}")
    if not os.path.exists(os.path.join(dir_name, "parameters")):
        os.symlink(f"{sim_dir}/parameters", os.path.join(dir_name, "parameters"))
    subprocess.run(f'{alf_dir}/Analysis/Max_SAC.out',
        check=True, cwd=dir_name, env={"OMP_NUM_THREADS": "1"})


def calc_k(alf_dir: Path, sim_dir: Path, results_dir: Path, always: bool = False) -> None:
    """Calculate k from k1, k2, k3 and write results to k.yaml.
    
    alf_dir: Path to ALF directory
    sim_dir: Path to simulation directory
    results_dir: Path where results should be written
    always: If True, always calculate k even if results already exist and are up to date
    """
    if not results_dir.exists():
        results_dir.mkdir(parents=True)
    if (results_dir / "k.yaml").exists() and (sim_dir / "data.h5").exists() and not always:
        if (sim_dir / "data.h5").stat().st_mtime < (results_dir / "k.yaml").stat().st_mtime:
            print(f"{results_dir / 'k.yaml'} is up to date, skipping k calculation.")
            return

    # Generate file Spintot
    subprocess.run(alf_dir / "Analysis" / "calc_k2_tau.out", check=True, cwd=sim_dir)
    os.rename(sim_dir / "Spintot", results_dir / "Spintot")

    # Generate replicas of Spintot.
    generate_tau_replicas(results_dir / "Spintot", 16)

    print(f"Running Max_SAC in {results_dir} for all Spintot replicas")
    dirs = list(results_dir.glob("Spintot_*"))

    with Pool() as pool:
        pool.map(partial(_run_max_sac, alf_dir, sim_dir), dirs)

    # Compute k2.
    k2_values = []
    for directory in dirs:
        k2 = calc_k2(directory)
        print(f"Calculated k2 for {directory}: {k2:.8E}")
        k2_values.append(k2)
    k2 = np.mean(k2_values)
    k2_err = np.std(k2_values)

    # Compute k1 and k3.
    k1, k1_err, k3, k3_err = calc_k1k3(sim_dir)

    # Combine k1, k2, k3 to get k.
    k_mean = k1 - (k2 - k3)
    k_err = np.sqrt(k1_err**2 + k2_err**2 + k3_err**2)

    params = f90nml.read(sim_dir / 'parameters')
    b_qmc = params["VAR_Kit_Heis"]["Hab"]
    theta_qmc = params["VAR_Kit_Heis"]["Htheta"]
    beta_qmc = params["VAR_Kit_Heis"]["Beta"]

    with open(sim_dir / "kitaev_K_abs", "r", encoding="utf-8") as f:
        kitaev_k_abs = float(f.readline().strip())

    mu_b = 0.05788
    with open(results_dir / "k.yaml", "w", encoding="utf-8") as f:
        yaml.dump({
            "B_QMC": abs(b_qmc),
            "theta_QMC": theta_qmc,
            "beta_QMC": beta_qmc,
            "B[T]": abs(b_qmc)*kitaev_k_abs/mu_b,
            "T[K]": kitaev_k_abs*11.6/beta_qmc,
            "k_av[mev]": float(k_mean*kitaev_k_abs),
            "k_err[mev]": float(k_err*kitaev_k_abs),
            "k1": float(k1),
            "k1_err": float(k1_err),
            "k2": float(k2),
            "k2_err": float(k2_err),
            "k3": float(k3),
            "k3_err": float(k3_err),
        }, f)

def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate the magnetotropic susceptibility k.")
    parser.add_argument("--alf-dir", type=Path, help="Path to ALF directory", required=True)
    parser.add_argument("--sim-dir", type=Path, default=Path.cwd(), help="Path to simulation directory")
    parser.add_argument("--results-dir", type=Path, default=Path.cwd() / "k", help="Path where analysis results should be written")
    parser.add_argument("--always", action="store_true", help="Always calculate k, even if results exist")
    args = parser.parse_args()
    if not args.results_dir.is_absolute():
        args.results_dir = args.sim_dir / args.results_dir
    calc_k(args.alf_dir, args.sim_dir, args.results_dir, always=args.always)

if __name__ == "__main__":
    main()
