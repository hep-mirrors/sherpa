#!/usr/bin/env python3
# Usage: ./plot-FF.py file1.csv file2.csv
# For example: ./plot-FF.py lead208_beam0_FF_q2.csv
# Generates line plots from CSV files for the Q2-dependent form factors.

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def plot_line(filepath):
    """Load CSV and create a 2D line plot."""
    try:
        # Load data
        data = np.loadtxt(filepath, delimiter=",", skiprows=1)
    except Exception as e:
        print(f"Error loading {filepath}: {e}", file=sys.stderr)
        return False

    if data.ndim == 1 or data.shape[1] < 2:
        print(
            f"Error: {filepath} must have at least 2 columns", file=sys.stderr
        )
        return False

    x, y = data[:, 0], data[:, 1]

    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot data
    ax.plot(x, y, linewidth=1.5)

    ax.set_xscale("log")
    ax.set_xlabel(r"$Q^2$ [GeV$^2$]")
    ax.set_ylabel(r"$F(Q^2)$")
    ax.grid(True, linestyle="--", alpha=0.7)

    # Save with same base name as input
    output_path = Path(filepath).with_suffix(".pdf")
    plt.tight_layout()
    plt.savefig(output_path, format="pdf")
    plt.close()

    print(f"Saved: {output_path}")
    return True


def main():
    if len(sys.argv) < 2:
        print("Usage: ./plot-FF.py file1.csv file2.csv ...", file=sys.stderr)
        sys.exit(1)

    success_count = 0
    for filepath in sys.argv[1:]:
        if plot_line(filepath):
            success_count += 1

    if success_count == 0:
        print("No plots generated successfully", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
