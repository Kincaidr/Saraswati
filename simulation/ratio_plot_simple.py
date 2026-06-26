#!/usr/bin/env python
# robert.kincaid@epfl.ch

import os

import matplotlib.pyplot as plt
import numpy as np


def compute_ratio(recovered_file, inj_sources=500, no_simulations=None, cap_counts=True):
    """
    Read recovered counts table and compute detection fraction per flux bin.

    File format (from radio_simulation.py):
      row 0: flux bin centres [mJy]
      rows 1+: recovered source counts per simulation

    ratio = sum(rec_counts) / (inj_sources * no_simulations)
    """
    data = np.loadtxt(recovered_file)
    bin_centres = data[0, :]

    if no_simulations is None:
        no_simulations = data.shape[0] - 1
    else:
        no_simulations = min(no_simulations, data.shape[0] - 1)

    rec_counts = np.sum(data[1 : 1 + no_simulations], axis=0)
    inj_total = inj_sources * no_simulations

    if cap_counts:
        rec_counts = np.minimum(rec_counts, inj_total)

    ratio = rec_counts / inj_total
    with np.errstate(divide="ignore", invalid="ignore"):
        error = ratio * np.sqrt(1.0 / inj_total + 1.0 / rec_counts)
    error = np.where(rec_counts > 0, error, np.nan)

    return bin_centres, ratio, error, rec_counts


def write_output_table(output_file, bin_centres, ratio, error):
    with open(output_file, "w") as file:
        for flux, rat, err in zip(bin_centres, ratio, error):
            file.write(f"{flux} {rat} {err}\n")
    print(f"Ratio table written to {output_file}")


def plot_ratio(bin_centres, ratio, error, output_plot):
    plt.errorbar(bin_centres, ratio, yerr=error, fmt="o", label="Detected fraction")
    plt.xlabel(r"Flux $S$ [mJy]", size=18)
    plt.ylabel("Detected fraction", size=18)
    plt.tick_params(axis="both", which="major", labelsize=15, length=5, width=1)
    plt.tick_params(axis="both", which="minor", labelsize=15, length=5, width=1)
    plt.axhline(y=1, color="green", linestyle="--")
    plt.xscale("log")
    plt.ylim(0, 1.05)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_plot)
    plt.show()
    plt.close()
    print(f"Plot saved to {output_plot}")


if __name__ == "__main__":
    name = "Zwcl2341"
    inj_sources = 500
    no_simulations = 20

    recovered_file = name + "_recovered_counts_table.txt"
    output_table = name + "_ratio_table.txt"
    plots_dir = "/home/kincaid/Desktop/Saraswati_codes/" + name + "/plots/"
    output_plot = os.path.join(plots_dir, name + "_ratio_plot.png")

    os.makedirs(plots_dir, exist_ok=True)

    bin_centres, ratio, error, rec_counts = compute_ratio(
        recovered_file,
        inj_sources=inj_sources,
        no_simulations=no_simulations,
    )

    print("Bin centres:", bin_centres)
    print("Summed recovered counts:", rec_counts)
    print("Ratio (rec / {}): {}".format(inj_sources * no_simulations, ratio))

    write_output_table(output_table, bin_centres, ratio, error)
    plot_ratio(bin_centres, ratio, error, output_plot)
