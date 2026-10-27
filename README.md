# CL20-HMX-NNP-MD

This repository contains computational input files, initial structures, a trained neural network potential, and analysis scripts associated with the manuscript:

**Thermal decomposition mechanism of a CL-20/HMX cocrystal revealed by neural network potential molecular dynamics simulations**

## Repository contents

The computational files are organized into four subdirectories within the uploaded project folder.

### init

Input files for initial structure preparation and reference calculations:

- `POSCAR`: initial atomic structure.
- `INCAR_rlx`: VASP settings for structural relaxation.
- `INCAR_md`: VASP settings for molecular dynamics.
- `param.json`: workflow parameters.
- `machine.json`: computing environment and job-submission configuration.

### run

Configuration files for the DP-GEN workflow and the trained potential:

- `INCAR_sp`: VASP settings for single-point reference calculations.
- `param.json`: workflow parameters.
- `machine.json`: computing environment and job-submission configuration.
- `frozen_model.pb`: trained neural network potential.

### md

LAMMPS simulation input files and initial structures for three systems:

- CL-20/HMX cocrystal.
- β-CL-20.
- γ-CL-20.

The files are organized into separate subdirectories for each system.

### analysis

- `compare_cl20_hbonds.py`: analyzes intermolecular C–H···O contacts and generates comparison plots with mean values and standard deviations.
- `species_evolution_cocrystal.py`: extracts time-dependent species counts for the CL-20/HMX cocrystal from a `dump.lammpstrj.species` file and writes the results to `data.txt`.

## Software requirements

The computational workflows use VASP, DP-GEN, DeePMD-kit, and LAMMPS. Species identification uses ReacNetGenerator.

The analysis scripts use Python. The contact-analysis script additionally requires NumPy and Matplotlib. These software packages must be installed separately.

## Usage notes

Computing environment settings in `machine.json` and file paths in the input files should be adapted to the local system.

The contact-analysis script requires trajectory and initial-structure files. Run `python compare_cl20_hbonds.py --help` to view its command-line options.

The species-extraction script should be run in a directory containing `dump.lammpstrj.species`. It writes `data.txt` in the working directory and replaces an existing file with that name.

Simulation inputs and analysis scripts are distinct from the trajectory files and species-identification outputs required for analysis.
