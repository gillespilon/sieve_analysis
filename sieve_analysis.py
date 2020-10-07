#! /usr/bin/env python3
"""
The purpose of this script is to perform a sieve analysis.

A sieve analysis is used to determine the particle size distribution of a
granular material. Material is passed through a series of progressively
smaller sieves. The mass of the material stopped by each sieve is determined
as a fraction of the whole mass.

[sieve_data.csv](https://drive.google.com/open?id=1QuhQmVAnxEakP879FnAt-GD0ElpuT49M)

This file has five columns:

- std sieve
- tyler sieve
- particle diameter
- sieve mass
- sieve soil mass

Clear the contents of the "sieve mass" and "sieve soil mass" columns.
Enter values of "sieve mass" and "sieve soil mass" for the desired sieves.
Save the file as a CSV with UTF-8 encoding.

To use the script, edit the constant for the raw material being sieved.
Execute the script.
"""

import pandas as pd
import numpy as np

colour1 = '#0077bb'
colour2 = '#33bbee'
# Enter the density of the material that was sieved.
density = 1.32
df = pd.read_csv('sieve_data.csv')
df['retained mass'] = df['sieve soil mass'] - df['sieve mass']
df = df.dropna(subset=['sieve mass'])
# rmt = retained mass total
rmt = df['retained mass'].sum()
df['retained pct'] = df['retained mass'] / rmt * 100
df['cumul retained pct'] = df['retained pct'].cumsum()
df['passing pct'] = 100 - df['cumul retained pct']
df.round(3)
df.to_csv('sieve_results.csv')
# Perform the geometric analysis
# gmps = geometric mean particle size
gmps = np.exp(((df['retained mass'] *\
       np.log(df['particle diameter'])).sum()) / rmt)
print('geometric mean particle size',
      np.exp(((df['retained mass'] *\
      np.log(df['particle diameter'])).sum())\
      /rmt).round(3),
      sep=" = ")
# gsd = geometric standard deviation
gsd = np.exp((((df['retained mass'] * (np.log(df['particle diameter']) -\
      (df['retained mass'] *\
      np.log(df['particle diameter'])).sum() / rmt)**2).sum()) / rmt)**0.5)
print('geometric standard deviation',
      gsd.round(3),
      sep=" = ")
surface_area = 6/density * np.exp(0.5 * (np.log(gsd)**2) - np.log(gmps / 10000))
print('surface area',
      surface_area.round(3),
      sep=" = ")
# nppg = number parts per g
nppg = 1/density * np.exp((4.5 * np.log(gsd)**2) - 3 * np.log(gmps / 10000))
print('number parts per g',
      nppg.round(3),
      sep=" = ")
# Create graphs
ax = df.plot(x='particle diameter', y='passing pct', logx=True, style='.-',\
             legend=False, color=colour1, grid=True)
ax.set_xlabel('Particle diameter (micron)')
ax.set_ylabel('Percent passing')
ax.set_title('Percent passing versus log particle diameter (%)')
ax.grid(True, which='minor', axis='x')
for spine in 'right', 'top':
    ax.spines[spine].set_visible(False)
ax.figure.savefig('sieve_percent_passing_vs_logdiameter.svg', format='svg')
ax.figure.savefig('sieve_percent_passing_vs_logdiameter.pdf', format='pdf')

ax = df.plot(x='particle diameter', y='passing pct', style='.-',\
             legend=False, color=colour2, grid=True)
ax.set_xlabel('Particle diameter (micron)')
ax.set_ylabel('Percent passing')
ax.set_title('Percent passing versus particle diameter (%)')
ax.grid(True, which='minor', axis='x')
for spine in 'right', 'top':
    ax.spines[spine].set_visible(False)
ax.figure.savefig('sieve_percent_passing_vs_diameter.svg', format='svg')
ax.figure.savefig('sieve_percent_passing_vs_diameter.pdf', format='pdf')
