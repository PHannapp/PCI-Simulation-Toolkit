# PCI Simulation Toolkit

## Description

The PCI Simulation Toolkit is an open-source Python framework developed to predict pressure–composition–temperature (PCT) diagrams of metal hydrides based solely on the atomic structure of the hydrogen-free metallic phase. It combines automated symmetry analysis, interstitial site detection, and configuration generation with ab initio DFT calculations to construct thermodynamic databases under para-equilibrium conditions. This toolkit supports the research presented in:

[10.XXXXXX](https://doi.org/10.XXXXXXX)

## Installation
~ Clone repo, installation of VASP, Enter POTCARS in dir POTCAR, install all dependencies of python libs

## Usage
Before using the toolkit, you need to install MongoDB and create two databases which names align with the settings in **MongoDB/connect.py**. 
Standard names are "System" for the systems database and "calculation" for the DFT database. 
The two json files provided (**DFTforCALPHAD.collection_calculation_pub.json** and **DFTforCALPHAD.System_pub.json**) can be readily imported in the MongoDB and used or extended.
The thermodynamic databases used in the publication are already stored in the systems database.

The program consists of five main Python files:

1. **_InsertSystem.py**: Insert a metall structure by importing a POSCAR file and decorate it with various elements. Here the interstitial sites will be determined.
2. **_preprocessing.py**:
    add_element() creates the POSCAR files for VASP calculations and stores them in the MongoDB
    vasp_input() creates all other necessary files (POTCAR, INCAR, KPOINTS) and stores them in the MongoDB
    process_output() processes the OUTCAR file from a finished DFT calculation
4. **_postprocessing.py**: calculates the formation enthalpies for every database entry
5. **_tdb_from_system.py**: creates a tdb ready to use for the para-equilibrium calculator stored in directory **CALPHAD_CALC**

## Citation
If you use this work, please cite it as:

@article{XX
}
