#! /usr/bin/env python3
"""
Perform a sieve analysis.

A sieve analysis is used to determine the particle size distribution of a
granular material. Material is passed through a series of progressively
smaller sieves. The mass of the material stopped by each sieve is determined
as a fraction of the whole mass.

[sieve_data.csv]
(https://drive.google.com/open?id=1QuhQmVAnxEakP879FnAt-GD0ElpuT49M)

This file has five columns:

- std sieve
- tyler sieve
- particle diameter
- sieve mass
- sieve soil mass

Clear the contents of the "sieve mass" and "sieve soil mass" columns.
Enter values of "sieve mass" and "sieve soil mass" for the desired sieves.
Save the file as a CSV with UTF-8 encoding.

To use the script, edit the density for the raw material being sieved.
Execute the script.
"""

from pathlib import Path
import time

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import datasense as ds


def main():
    start_time = time.perf_counter()
    path_logdiameter = "sieve_percent_passing_vs_logdiameter.svg"
    path_diameter = "sieve_percent_passing_vs_diameter.svg"
    path_file_out = Path("sieve_results.csv")
    path_file_in = Path("sieve_data.csv")
    output_url = "sieve_analysis.html"
    header_title = "Sieve Analysis"
    header_id = "sieve_analysis"
    figsize = (8, 6)
    grid_alpha = 0.2
    density = 1.32
    original_stdout = ds.html_begin(
        output_url=output_url, header_title=header_title, header_id=header_id
    )
    ds.script_summary(script_path=Path(__file__), action="started at")
    ds.style_graph()
    df = pd.read_csv(path_file_in)
    df["retained mass"] = df["sieve soil mass"] - df["sieve mass"]
    df = df.dropna(axis=0, how="any", subset=["sieve mass"])
    retained_mass_total = df["retained mass"].sum()
    df["retained percentage"] = df["retained mass"] / retained_mass_total * 100
    df["cumulative retained percentage"] = df["retained percentage"].cumsum()
    df["passing percentage"] = 100 - df["cumulative retained percentage"]
    df.round(3).to_csv(path_file_out)
    geometric_mean_particle_size = np.exp(
        ((df["retained mass"] * np.log(df["particle diameter"])).sum()) /
        retained_mass_total
    )
    print(f"geometric mean particle size {geometric_mean_particle_size:10.3f}")
    geometric_standard_deviation = np.exp(
        (((df["retained mass"] * (np.log(df["particle diameter"]) -
         (df["retained mass"] * np.log(df["particle diameter"])).sum() /
         retained_mass_total)**2).sum()) / retained_mass_total)**0.5
    )
    print(f"geometric standard deviation {geometric_standard_deviation:10.3f}")
    surface_area = 6 / density * np.exp(
        0.5 * (np.log(geometric_standard_deviation)**2) -
        np.log(geometric_mean_particle_size / 10000))
    print(f"surface area                 {surface_area:10.3f}")
    number_parts_per_g = 1 / density * \
        np.exp(
            (4.5 * np.log(geometric_standard_deviation)**2) -
            3 * np.log(geometric_mean_particle_size / 10000)
        )
    print(f"number parts per g           {number_parts_per_g:10.3f}")
    print()
    fig, ax = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=figsize
    )
    fig.suptitle(t="Sieve Analysis")
    ax.semilogx(df["particle diameter"], df["passing percentage"])
    ax.grid(visible=True, which="both", axis="both", alpha=grid_alpha)
    ax.set_title(label="Percent passing versus log particle diameter")
    ax.set_xlabel(xlabel="Particle diameter (micron)")
    ax.set_ylabel(ylabel="Passing (%)")
    ds.despine(ax=ax)
    fig.savefig(fname=path_logdiameter, format="svg")
    ds.html_figure(file_name=path_logdiameter)
    fig, ax = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=figsize
    )
    fig.suptitle(t="Sieve Analysis")
    ax.plot(df["particle diameter"], df["passing percentage"])
    ax.grid(visible=True, which="both", axis="both", alpha=grid_alpha)
    ax.set_title(label="Percent passing versus particle diameter")
    ax.set_xlabel(xlabel="Particle diameter (micron)")
    ax.set_ylabel(ylabel="Passing (%)")
    ds.despine(ax=ax)
    fig.savefig(fname=path_diameter, format="svg")
    ds.html_figure(file_name=path_diameter)
    stop_time = time.perf_counter()
    ds.script_summary(script_path=Path(__file__), action="finished at")
    ds.report_summary(start_time=start_time, stop_time=stop_time)
    ds.html_end(original_stdout=original_stdout, output_url=output_url)


if __name__ == "__main__":
    main()
