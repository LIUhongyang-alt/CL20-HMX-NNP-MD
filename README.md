# CL20-HMX-NNP-MD

This repository contains computational input files, initial structures, and a trained neural network potential associated with the manuscript:

**Thermal decomposition mechanism of a CL-20/HMX cocrystal revealed by neural network potential molecular dynamics simulations**

## Repository contents

The computational files are organized into three subdirectories within the `cl20-hmx` folder.

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

## Software requirements

The computational workflows use VASP, DP-GEN, DeePMD-kit, and LAMMPS. These software packages must be installed separately.

## Usage notes

Computing environment settings in `machine.json` and file paths in the input files should be adapted to the local system.

The trained neural network potential is provided as `frozen_model.pb` for use with a compatible DeePMD-kit and LAMMPS installation.
