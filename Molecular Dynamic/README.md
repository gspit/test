# Polymer LAMMPS Molecular Dynamics Simulation Results

## Project Overview

This directory contains LAMMPS molecular dynamics simulation results for 6 polymer candidate systems, used to calculate melting temperature (Tm).

---

## Folder Overview

|  Folder  |  Status | | Availability |
|--------  | --------|-|--------------|
| ANNTm1 | ✅ Complete | ✅ Recommended |
| ANNTm2 | ✅ Complete | ✅ Recommended |
| ANNTm3 | ✅ Complete | ✅ Recommended |
| DESCTm1 | ✅ Complete |✅ Recommended |
| DESCTm2 | ✅ Complete |✅ Recommended |
| DESCTm3 | ✅ Complete |✅ Recommended |

## File Description

### Core Simulation Files

| File | Description |
|------|-------------|
| **in2.txt** | LAMMPS input script — main control file with detailed comments, defines simulation parameters (ensemble, temperature ramp, pressure control, etc.) |
| **lammps.data** | Initial structure file — contains atomic coordinates, bond/angle/dihedral topology, force field parameters |
| **log.lammps** | Simulation log — complete output including energy, temperature, pressure, volume at each timestep |
| **H_T.dat** | Enthalpy-Temperature data — collected during heating scan, used for Tg/Tm calculation |

### Analysis Files

| File | Description |
|------|-------------|
| **fit.py** | Python fitting script — analyzes H_T.dat to extract Tg/Tm via statistical mechanical methods |
| **test.py** | Alternative analysis script |
| **Tm.png** | Melting temperature fitting plot |
| **fit.png** | Fitting curve visualization |
| **test.png** | Test/validation plot |

### Root Directory

| File | Description |
|------|-------------|
| **lammps.job** | SLURM job submission script — used for HPC cluster submission, not needed for local runs |

## Required Files for Thesis/Publication

1. **Input**: `in2.txt`
2. **Structure**: `lammps.data`
3. **Results**: `log.lammps`, `H_T.dat`
4. **Figures**: `Tm.png`
5. **Analysis**: `fit.py` (for reproducing fitting results)

---

## Acknowledgments

We gratefully acknowledge the computing resources provided by:

- **Hebei University** 
- **Nanjing University of Aeronautics and Astronautics** 

The simulations were performed on the supercomputing clusters of both institutions.
