from Equations import (
    Structure,
    GLOBAL_Structure,
    gibbs_minimizer,
)  # Import custom classes and functions for thermodynamic calculations
import Calculate_PCI  # Module for calculating Pressure-Composition Isotherms
import numpy as np  # For numerical operations
import matplotlib.pyplot as plt  # For plotting
import random
import time

database_modules = ["database"] 
name = ['alpha']
solids = [[["TI", "ZR"], ["MN"]]]  # Define solids involved in the structure
multiplicities = [[1,2,1,1,2]]  # Define the multiplicity of sites in the structure
x_La = 1
x_Cu = 0
site_fractions = [[[x_La, 1 - x_La], [x_Cu, 1 - x_Cu]]]  # Initial site fractions for hydrogen
color = ["blue"]
MIN_Y = 0
MAX_Y = 24
x_unit = "HM"

vars = [[]]

dir = '/home/pet69516/Dokumente/DFTforCALPHAD/CALPHAD_CALC/'
T = 296  # Set temperature (in Kelvin)
t0 = time.time()
for i in range(1):
    structure_list = []

    for i in range(len(database_modules)):
        module_name = database_modules[i]
        multiplicity = multiplicities[i]
        solid = solids[i]
        site_fraction = site_fractions[i]
        # Create a structure instance
        Initial_Structure = Structure(solid, multiplicity, site_fraction, name[i], module_name)
        structure_list.append(Initial_Structure)
        # Perform Gibbs energy minimization to find equilibrium compositions --------------------------
        vars_part = vars[i]
        gibbs_minimizer(T, Initial_Structure, vars_part) 

    for i, structure in enumerate(structure_list):
        Calculate_PCI.plateau(structure, T) 

        # if no plateaus were found -> search for second order plateaus (from minima or ends)
        if len(structure.plateau_pressures) == 0:
            Calculate_PCI.plateau_from_minima(structure, T)

    if len(structure_list) > 0:
        # Create a structure instance
        global_structure = GLOBAL_Structure(solids, multiplicities, site_fraction, name, database_modules)  

        Calculate_PCI.plateau_between_phases(global_structure, structure_list, T)   

        result, slope, intercept = Calculate_PCI.calculate_x_from_p(global_structure, 1e5, T)
        
    print(f"Time used for 1000 equilibria: {time.time() - t0} s.")

    for plateau in range(global_structure.n_plateau):
        print(f"{name[global_structure.phase_id[plateau][0]]} is stable until {global_structure.x_H[global_structure.CTs[plateau][0]]:.3f}.")
        print(f"{name[global_structure.phase_id[plateau][1]]} is stable from {global_structure.x_H[global_structure.CTs[plateau][1]]:.3f}.")
        print(f"The plateau pressure of the {plateau + 1}. plateau is at {global_structure.plateau_pressures[plateau] / 1e5:.2f} bar.")
    
    plt.figure()
    plt.plot(global_structure.x_H, global_structure.gm)

    for m, CT in enumerate(global_structure.CTs):
        plt.plot([global_structure.x_H[CT[0]], global_structure.x_H[CT[1]]], [global_structure.gm[CT[0]], global_structure.gm[CT[1]]], 'o')
        CT_x = np.linspace(global_structure.x_H[CT[0]], global_structure.x_H[CT[1]], 100)
        CT_gm = global_structure.plateau_slopes[m] * CT_x + global_structure.plateau_intercepts[m]  
        plt.plot(CT_x, CT_gm)

    pl_root = []
    p_plot = np.array([])
    if global_structure.CTs[0][0] != global_structure.CTs[0][1]:
        for p, plateau in enumerate(global_structure.plateau_pressures):
                pl_root.append(np.log10(
                    plateau
                ))  # Calculate the logarithmic base for plateau pressures
                # Adjust pressure points for detailed plotting around plateau pressures
                p_plot = np.append(p_plot, plateau + plateau / 200000000)
                p_plot = np.append(p_plot, plateau - plateau / 200000000)
        for plateau in range(len(global_structure.plateau_pressures) - 1):
            p_plot = np.append(p_plot, np.logspace(
                pl_root[plateau], pl_root[plateau + 1], 200
            ))  # Generate pressures for plotting
        min = pl_root[0] - 3
        if min > MIN_Y:
            min = MIN_Y
        p_plot = np.append(p_plot, np.logspace(
            min, pl_root[0], 200
        ))  # Generate pressures for plotting
        max = pl_root[-1] + 3
        if max < MAX_Y:
            max = MAX_Y
        p_plot = np.append(p_plot, np.logspace(
            pl_root[-1], max, 200
        ))  # Generate pressures for plotting
    else:
        p_plot = np.logspace(MIN_Y,MAX_Y,200)
    p_plot = np.sort(p_plot)  # Sort the pressures for plotting
    x_plot = []  # Initialize list for hydrogen mole fractions
    # Calculate hydrogen mole fraction for each pressure point
    for p in p_plot:
        result, slope, intercept = Calculate_PCI.calculate_x_from_p(global_structure, p, T)
        x_plot.append(global_structure.x_H[result])
    plt.savefig(f'{dir}/plot_gibbs.png')  
    plt.figure(3)
    if x_unit == "HM":
        x_plot = [x / (1-x) for x in x_plot]

    plt.semilogy(
        x_plot, p_plot, 
        label=f"{global_structure.name}", color=color[i]
    )  # Create a semilogarithmic plot
    plt.ylabel("Hydrogen Pressure, $\t{P}$ / Pa")
    plt.xlabel("Hydrogen mole fraction, $\t{x_{H}}$ / ")
    plt.tight_layout()  # Adjust
    plt.legend()
    plt.grid()
    plt.savefig(f'{dir}/plot_PCI.png')
    #plt.show()