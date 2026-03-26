# Prediction of Polymer Melting Temperature Using Deep Learning, Descriptor Models, and Molecular Dynamics Validation

## Overview

This repository contains the complete implementation, data, and validation results for predicting polymer melting temperature (Tm) using these approaches:


## Repository Structure

```
.
├── data/                          # Training and test datasets
│   ├── train_data.csv             # 3,660 polymers with Tm
│   ├── Regression_Tm.xlsx         # Regression analysis dataset
│   ├── PI1M50K_Prediction_Results.xlsx  # Model predictions
│   ├── TableS1S2S3.xlsx           # Supplementary tables
│   └── Candidates_from_PI1M.csv   # 54,518 Candidate polymers from PI1M database
│   ├── stepwise_regression.py     # Stepwise regression implementation (NumPy + SciPy)
│   └── README_Stepwise_Regression.txt   # Regression model documentation
│
├── molecular_dynamics/            # MD simulation data and scripts
│   ├── lammps.job                 # LAMMPS job submission script
│   ├── ANNTm1/, ANNTm2/, ANNTm3/  # Deep learning model validation (3 replicates)
│   ├── DESCTm1/, DESCTm2/, DESCTm3/    # Descriptor model validation (3 replicates)
│   └── README.md                  # MD simulation documentation

```


## Models

### 1. Stepwise Regression (stepwise_regression.py)

Pure Python implementation using NumPy and SciPy (no external ML libraries required).

**Usage**:
```bash
python stepwise_regression.py
```

### 2. Deep Learning Models

This study adopts the **Uni-QSAR platform** (https://www.bohrium.com/apps/qsar-web-new) — a QSAR modeling tool provided by Bohrium Cloud — for predicting polymer melting temperature:
- Input: Polymer SMILES strings as the only explicit input (no manual molecular feature extraction required)
- Modeling Framework: An end-to-end automated multi-modal learning framework that can automatically derive and learn multi-level molecular structure information from SMILES representations (not limited to one-dimensional sequence features)

## Molecular Dynamics Validation

### LAMMPS Simulations

Each model (ANN and DESC) is validated on 3 polymer systems using LAMMPS:

**Directory Structure**:
```
molecular_dynamics/
├── ANNTm1/, ANNTm2/, ANNTm3/    # Deep learning validation
├── DESCTm1/, DESCTm2/, DESCTm3/ # Descriptor model validation
└── Each contains:
    ├── fit.py                   # Python fitting script
    ├── in2.txt                  # LAMMPS input file
    ├── lammps.data              # Molecular structure data
    ├── H_T.dat                  # Enthalpy vs. Temperature data
    ├── log.lammps               # LAMMPS simulation log
    └── Tm.png                   # Melting point visualization
```

### Python Dependencies
```bash
pip install numpy scipy pandas
```

### LAMMPS (for MD simulations)
- Download: https://www.lammps.org/
- Installation: Follow official documentation for your OS

### Data Files
All data files are included in the `data/` directory. No external downloads required.



## Citation

If you use this code or data in your research, please cite:

```bibtex
@article{Guo2026,
  title={Prediction of Polymer Melting Temperature Using Deep Learning, Descriptor Models, and Molecular Dynamics Validation},
  author={Guo, Haiqian and Wang, Yaling and Zan, Wei and Liu, Jiejie and Ba, Xinwu},
  journal={Computational Materials Science},
  year={2026},
  volume={},
  pages={},
  doi={}
}
```

## License

This project is licensed under the MIT License — see the LICENSE file for details.
