# PCI Simulation Toolkit

## Description

The PCI Simulation Toolkit is an open-source Python framework developed to predict pressure–composition–temperature (PCT) diagrams of metal hydrides based solely on the atomic structure of the hydrogen-free metallic phase. It combines automated symmetry analysis, interstitial site detection, and configuration generation with ab initio DFT calculations to construct thermodynamic databases under para-equilibrium conditions. This toolkit supports the research presented in:

[10.XXXXXX](https://doi.org/10.XXXXXXX)

## Installation
~ Clone repo, installation of VASP, Enter POTCARS in dir POTCAR, install all dependencies of python libs

## Usage
The program consists of five main Python files:

1. **Calculate_PCI.py**: Calculates the plateau pressures and pressure values corresponding to mole fractions in single-phase regions.
2. **CONSTANTS.py**: Contains constants used for calculating the common tangent.
3. **Equations.py**: Minimizes the Gibbs energy curve of the solid phase by finding the global minimum for discrete compositions of xH.
4. **database.py**: Contains the thermodynamic variables to describe different phases.
5. **Execute.py**: Provides example usage of the files to generate PCI diagrams.

## Running the Code
To run the program, the file `Execute.py` gives example usage. The script will import the necessary files and prompt the user to input the following variables:

- **solids**: A list of n lists representing the composition of each sublattice.
- **multiplicity**: A list of integers defining the multiplicity of the sublattices.
- **site_fractions**: A list of n lists representing the site fractions of the metal elements in each sublattice.

## Outputs
Running Execute.py will generate the PCI diagrams as described in Sections S1-S3 of the supplementary material.

## Citation
If you use this work, please cite it as:

@article{XX
}
