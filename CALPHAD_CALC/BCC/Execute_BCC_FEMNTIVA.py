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
import os

dir = r"/home/pet69516/Dokumente/DFTforCALPHAD/CALPHAD_CALC/BCC/"

# Choose a colormap
cmap = plt.cm.get_cmap('Set1', 10)
file_names_BCC = [
#"database_BCC_10_nothing",
"database_BCC_1",
#"database_BCC_10_both",
]
file_names_FCC = [
#"database_FCC_both",
"database_FCC_2",
#"database_FCC_both",
]
names_BCC = [
"BCC",
#"PBE_Entropy",
]
names_FCC = [
"FCC",
#"PBE_Entropy",
]
multiplicitiess_BCC = [
[1,1]
]
multiplicitiess_FCC = [
[1,2]
]
plot_descriptor = "BCC"

exp_colors = ["black"]
exp=[]
exp2= [
[
[0.06819361614795855, 128361.19008195761],
[0.07740782770359994, 316227.7660168379],
[0.09855554482458889, 513858.2972943202],
[0.11462080788082185, 598561.0238654562],
[0.2529184875775648, 659594.6257418732],
[0.4603378443634257, 673462.2438755755],
[0.5455142212135113, 692406.9936029648],
[0.6537306540792429, 692406.9936029648],
[0.7950177463360457, 706964.4739183798],
[0.8972342274911271, 747298.2683491874],
[1.05187483400536, 3436723.3870942583],
[1.0740758625685105, 7069644.7391837975],
[1.0913589226646065, 22201664.839136124],
[1.0871562401912258, 9074680.12178283],
[1.0618888379167977, 3342692.181383907],
[1.0316355602771812, 1375796.5525531906],
[0.9703956104981049, 800966.6945806034],
[0.9051702803196755, 517434.5473123614],
[0.7318175507641791, 503277.16864677577],
[0.590525931380834, 482764.06972582516],
[0.4532544003669992, 489507.1459699249],
[0.25785909167733057, 479427.44479029614],
[0.09149926358741578, 423161.37934508896],
[0.06830981239587605, 218958.31676150055],
]
]

exp_names = ["Experiment"]
for num in range(len(file_names_BCC)):
    file_name_BCC = file_names_BCC[num]
    file_name_FCC = file_names_FCC[num]

    database_modules = [file_name_BCC, file_name_FCC] 
    name = [names_BCC[num], names_FCC[num]]
    solids =[[["FE", "MN", "TI", "VA"]], [["FE", "MN", "TI", "VA"]]]  # Define solids involved in the structure
    multiplicities = [multiplicitiess_BCC[num], multiplicitiess_FCC[num]]  # Define the multiplicity of sites in the structure
    x_FE = 0.15
    x_MN = 0.05
    x_TI = 0.4   
    site_fractions = [[[x_FE, x_MN, x_TI, 1-x_FE-x_MN-x_TI]], [[x_FE, x_MN, x_TI, 1-x_FE-x_MN-x_TI]]]  # Initial site fractions for hydrogen
    print(1-x_FE-x_MN-x_TI)
    color = ["blue"]
    FORCE = True
    MIN_Y = -10
    MAX_Y = 10
    MIN_X = 0
    MAX_X = 2
    x_unit = "HM"

    vars = [[], []]

    T = 353  # Set temperature (in Kelvin)
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
            print(f"The plateau pressure of the {plateau + 1}. plateau is at {global_structure.plateau_pressures[plateau] / 1e5:.2e} bar.")

        plt.figure(2)
        plt.plot(global_structure.x_H, global_structure.gm, label=f"{global_structure.name[0]}", color=cmap(num))

        for m, CT in enumerate(global_structure.CTs):
            plt.plot([global_structure.x_H[CT[0]], global_structure.x_H[CT[1]]], [global_structure.gm[CT[0]], global_structure.gm[CT[1]]], 'o', color=cmap(num))
            CT_x = np.linspace(global_structure.x_H[CT[0]], global_structure.x_H[CT[1]], 100)
            CT_gm = global_structure.plateau_slopes[m] * CT_x + global_structure.plateau_intercepts[m]  
            plt.plot(CT_x, CT_gm, color=cmap(num))

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

        plt.figure(3)
        if x_unit == "HM":
            x_plot = [x / (1-x) for x in x_plot]
            plt.xlabel("Hydrogen atoms / Metal atom, $\t{HM}$ / ")
        else:
            plt.xlabel("Hydrogen mole fraction, $\t{x_{H}}$ / ")

        plt.semilogy(
            x_plot, p_plot, 
            label=f"{global_structure.name[0]}", color=cmap(num)
        )  # Create a semilogarithmic plot

if FORCE:
    plt.figure(3)
    plt.ylim(np.power(10.0,MIN_Y), np.power(10.0,MAX_Y))
    plt.xlim(MIN_X, MAX_X)
plt.figure(2)
plt.xlabel("Hydrogen mole fraction, $\t{x_{H}}$ / ")
plt.ylabel("Molar Gibbs Energy, $\t{G}$ / J / mol")
plt.tight_layout()  # Adjust
plt.legend()
plt.grid()
plt.savefig(f'{dir}/plot_gibbs_{plot_descriptor}.png')  
plt.close()
plt.figure(3)
for i in range(len(exp)):
    plt.plot([line[0] for line in exp[i]], [line[1] for line in exp[i]], 'o', color = exp_colors[i], linestyle='-', label=exp_names[i])
plt.ylabel("Hydrogen Pressure, $\t{P}$ / Pa")
plt.tight_layout()  # Adjust
plt.legend(loc='lower right', ncol=2)
plt.grid()
plt.savefig(f'{dir}/plot_PCI_{plot_descriptor}.png')
#plt.show()
plt.close()