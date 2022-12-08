#! /usr/bin/env python3
"""
Perform a sieve analysis.

A sieve analysis is used to determine the particle size distribution of a
granular material. Material is passed through a series of progressively
smaller sieves. The mass of the material stopped by each sieve is determined
as a fraction of the whole mass.
"""

from pathlib import Path
import time

import matplotlib.pyplot as plt
import datasense as ds
import pandas as pd
import numpy as np


def main():
    start_time = time.perf_counter()
    cumulative_retained_percentage = "cumulative retained percentage"
    path_logdiameter = "sieve_percent_passing_vs_logdiameter.svg"
    path_diameter = "sieve_percent_passing_vs_diameter.svg"
    retained_percentage = "retained percentage"
    passing_percentage = "passing percentage"
    path_file_out = Path("sieve_results.csv")
    particle_diameter = "particle diameter"
    path_file_in = Path("sieve_data.csv")
    sieve_soil_mass = "sieve soil mass"
    output_url = "sieve_analysis.html"
    header_title = "Sieve Analysis"
    retained_mass = "retained mass"
    header_id = "sieve_analysis"
    sieve_mass = "sieve mass"
    figsize = (8, 6)
    grid_alpha = 0.2
    density = 1.32
    original_stdout = ds.html_begin(
        output_url=output_url,
        header_title=header_title,
        header_id=header_id
    )
    ds.script_summary(
        script_path=Path(__file__),
        action="started at"
    )
    ds.style_graph()
    df = pd.DataFrame(
        data={
            particle_diameter: [4760, 2000, 841, 420, 250, 74, 10],
            sieve_mass: [1, 1, 208, 191.4, 72.3, 23.1, 0.9],
            sieve_soil_mass: [1.6, 1.6, 416, 382.8, 144.6, 46.2, 1.8]
        }
    )
    # df = ds.read_file(file_name=path_file_in)
    df[retained_mass] = df[sieve_soil_mass] - df[sieve_mass]
    # drop empty rows
    df = df[df[[sieve_mass]].notna().all(axis="columns")]
    retained_mass_total = df[retained_mass].sum()
    df[retained_percentage] = df[retained_mass] / retained_mass_total * 100
    df[cumulative_retained_percentage] = df[retained_percentage].cumsum()
    df[passing_percentage] = 100 - df[cumulative_retained_percentage]
    # ds.save_file(
    #     df=df.round(3),
    #     file_name=path_file_out
    # )
    geometric_mean_particle_size = np.exp(
        ((df[retained_mass] * np.log(df[particle_diameter])).sum()) /
        retained_mass_total
    )
    print(f"geometric mean particle size {geometric_mean_particle_size:10.3f}")
    geometric_standard_deviation = np.exp(
        (((df[retained_mass] * (np.log(df[particle_diameter]) -
         (df[retained_mass] * np.log(df[particle_diameter])).sum() /
         retained_mass_total)**2).sum()) / retained_mass_total) ** 0.5
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
    ax.semilogx(
        df[particle_diameter],
        df[passing_percentage]
    )
    ax.grid(
        visible=True,
        which="both",
        axis="both",
        alpha=grid_alpha
    )
    ax.set_title(label="Percent passing versus log particle diameter")
    ax.set_xlabel(xlabel="Particle diameter (micron)")
    ax.set_ylabel(ylabel="Passing (%)")
    ds.despine(ax=ax)
    fig.savefig(
        fname=path_logdiameter,
        format="svg"
    )
    ds.html_figure(file_name=path_logdiameter)
    fig, ax = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=figsize
    )
    fig.suptitle(t="Sieve Analysis")
    ax.plot(
        df[particle_diameter],
        df[passing_percentage]
    )
    ax.grid(
        visible=True,
        which="both",
        axis="both",
        alpha=grid_alpha
    )
    ax.set_title(label="Percent passing versus particle diameter")
    ax.set_xlabel(xlabel="Particle diameter (micron)")
    ax.set_ylabel(ylabel="Passing (%)")
    ds.despine(ax=ax)
    fig.savefig(
        fname=path_diameter,
        format="svg"
    )
    ds.html_figure(file_name=path_diameter)
    stop_time = time.perf_counter()
    ds.script_summary(
        script_path=Path(__file__),
        action="finished at"
    )
    ds.report_summary(
        start_time=start_time,
        stop_time=stop_time
    )
    ds.html_end(
        original_stdout=original_stdout,
        output_url=output_url
    )


if __name__ == "__main__":
    main()
