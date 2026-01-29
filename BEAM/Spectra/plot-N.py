#!/usr/bin/env python3
# Usage: ./plot-N.py file1.csv file2.csv ...
# For example: ./plot-N.py lead208_beam0_N_x_b.csv
# Generates contour plots with logarithmic x, b, and N axes from CSV files.
# Each CSV must have exactly 3 columns: x, b, N. The header is ignored.

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


def plot_csv(filepath):
    """Load CSV and create logarithmic contour plot."""
    try:
        data = np.loadtxt(filepath, delimiter=",", skiprows=1)
    except Exception as e:
        print(f"Error loading {filepath}: {e}", file=sys.stderr)
        return False

    if data.shape[1] != 3:
        print(
            f"Error: {filepath} must have exactly 3 columns, found {data.shape[1]}",
            file=sys.stderr,
        )
        return False

    x, y, z = data[:, 0], data[:, 1], data[:, 2]

    # Handle 1D edge case where y (column 2) is all zeros
    if np.allclose(y, 0):
        fig, ax = plt.subplots(figsize=(8, 6))

        # Sort by x to ensure a clean line plot
        sort_idx = np.argsort(x)
        ax.plot(x[sort_idx], z[sort_idx])

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("x")
        ax.set_ylabel(r"$N(x)$")
        ax.grid(True, linestyle="--", alpha=0.7)

        output_path = Path(filepath).with_suffix(".pdf")
        plt.tight_layout()
        plt.savefig(output_path, format="pdf")
        plt.close()

        print(f"Saved (1D mode): {output_path}")
        return True

    # Check for non-positive values (incompatible with log scale)
    if np.any(x <= 0) or np.any(y <= 0) or np.any(z <= 0):
        print(
            f"Warning: {filepath} contains non-positive values; log scale may fail",
            file=sys.stderr,
        )

    # Determine if data is on regular grid
    x_unique = np.unique(x)
    y_unique = np.unique(y)

    fig, ax = plt.subplots(figsize=(8, 6))

    # Filter strictly positive, finite values to determine color limits
    z_valid = z[(z > 0) & np.isfinite(z)]

    if z_valid.size == 0:
        print(
            f"Error: {filepath} has no valid positive data for log scale.",
            file=sys.stderr,
        )
        return False

    z_min, z_max = z_valid.min(), z_valid.max()

    # Handle edge case where data is constant
    if z_min == z_max:
        z_min, z_max = z_min * 0.9, z_max * 1.1

    # Generate log-spaced levels based on valid range
    levels = np.logspace(np.log10(z_min), np.log10(z_max), 18)

    if len(x_unique) * len(y_unique) == len(x):
        # Regular grid: reshape and use contourf
        nx, ny = len(x_unique), len(y_unique)
        X = x.reshape(ny, nx)
        Y = y.reshape(ny, nx)
        Z = z.reshape(ny, nx)
        contour = ax.contourf(X, Y, Z, levels=levels, norm=LogNorm())
    else:
        # Scattered points: use tricontourf
        contour = ax.tricontourf(x, y, z, levels=levels, norm=LogNorm())

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("x")
    ax.set_ylabel(r"$b$ [1/GeV]")

    # Generates 6 evenly spaced ticks in log space between min and max
    cbar_ticks = np.geomspace(z_min, z_max, 6)

    cbar = plt.colorbar(contour, ax=ax, ticks=cbar_ticks, format="%.1e")
    cbar.set_label(r"$N(x, b)$ [GeV]")

    # Save with same base name as input
    output_path = Path(filepath).with_suffix(".pdf")
    plt.tight_layout()
    plt.savefig(output_path, format="pdf")
    plt.close()

    print(f"Saved: {output_path}")
    return True


def main():
    if len(sys.argv) < 2:
        print("Usage: ./plot-N.py file1.csv file2.csv ...", file=sys.stderr)
        sys.exit(1)

    success_count = 0
    for filepath in sys.argv[1:]:
        if plot_csv(filepath):
            success_count += 1

    if success_count == 0:
        print("No plots generated successfully", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
