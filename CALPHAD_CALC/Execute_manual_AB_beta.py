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

# Choose a colormap
cmap = plt.cm.get_cmap('Set1', 10)
file_names = [
"database_manual_AB_beta",
#"database_TiMn2_no",
#"database_TiMn2_metal",
]
names = [
"DFT calculated",
#"No interactions",
#"Metal interaction",
]
multiplicitiess = [
 [[8, 8, 31]]  
]
plot_descriptor = "AB_beta_manual"

exp_colors = ["black", "green"]
exp= [
[
[0.016826923076923073, 4.761872663008537],
[0.04086538461538458, 8.504124891946162],
[0.06490384615384615, 11.17244301313065],
[0.09615384615384609, 13.250193550567484],
[0.18749999999999994, 15.187331811625278],
[0.3870192307692308, 15.984671064343198],
[0.5865384615384617, 15.984671064343198],
[0.6730769230769229, 18.01173528334133],
[0.7427884615384615, 21.914962948968714],
[0.8052884615384615, 29.53727118810853],
[0.8677884615384615, 37.503696379226895],
[0.939903846153846, 51.85788673392234],
[0.9783653846153846, 68.12920690579611],
[0.9907692307692306, 470.61872663008535],
[0.9807692307692306, 47.61872663008535],
[0.939903846153846, 23.065504755640514],
[0.8389423076923077, 18.636758192040272],
[0.75, 15.580901878706047],
[0.6057692307692306, 11.364636663857242],
[0.5192307692307692, 6.812920690579611],
[0.4543269230769231, 7.35642254459641],
[0.2836538461538461, 7.049334974921203],
[0.1490384615384615, 6.989473207273483],
[0.04567307692307687, 4.447829676127633],
[0.05288461538461531, 2.5769803745148776],
[0.026442307692307654, 1.8165997883753278],
]
]

for e in exp[0]:
    print(e)
    e[1] = e[1]*1e5
x_Tis = [1]
Ts = [273]
exp_names = ["FeTi - 273 K"]#, "LaNi4.5Al0.5", "LaNi4Al"]
for num in range(len(file_names)):
    for xps in range(len(x_Tis)):
        file_name = file_names[num]

        database_modules = [file_name] 
        name = [names[num]]
        solids = [[["FE"], ["TI"]]]  # Define solids involved in the structure
        multiplicities = multiplicitiess[num]  # Define the multiplicity of sites in the structure
        x_Ti = x_Tis[xps]
        site_fractions = [[[1], [1]]]  # Initial site fractions for hydrogen
        T = Ts[xps]  # Set temperature (in Kelvin)
        color = ["blue"]
        FORCE = True
        MIN_Y = 4
        MAX_Y = 9.7
        x_unit = "HM"

        vars = [[]]

        dir = '/home/pet69516/Dokumente/DFTforCALPHAD_home/CALPHAD_CALC/'
        t0 = time.time()
        for i in range(1):
            structure_list = []

            for i in range(len(database_modules)):
                module_name = database_modules[i]
                multiplicity = multiplicities[i]
                solid = solids[i]
                site_fraction = site_fractions[i]
                # Create a structure instance
                print(module_name)
                Initial_Structure = Structure(solid, multiplicity, site_fraction, name[i], module_name)
                structure_list.append(Initial_Structure)
                # Perform Gibbs energy minimization to find equilibrium compositions --------------------------
                vars_part = vars[i]
                gibbs_minimizer(T, Initial_Structure, vars_part, steps=100) 

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

            plt.figure(2)
            plt.plot(global_structure.x_H, global_structure.gm, label=f"{global_structure.name[0]}", color=cmap(num+xps))

            for m, CT in enumerate(global_structure.CTs):
                #plt.plot([global_structure.x_H[CT[0]], global_structure.x_H[CT[1]]], [global_structure.gm[CT[0]], global_structure.gm[CT[1]]], 'o', color=cmap(num))
                CT_x = np.linspace(global_structure.x_H[CT[0]], global_structure.x_H[CT[1]], 100)
                CT_gm = global_structure.plateau_slopes[m] * CT_x + global_structure.plateau_intercepts[m]  
                #plt.plot(CT_x, CT_gm, color=cmap(num))

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
                label=f"{global_structure.name[0]}", color=cmap(num+xps)
            )  # Create a semilogarithmic plot

if FORCE:
    plt.figure(3)
    plt.ylim(np.power(10.0,MIN_Y), np.power(10.0,MAX_Y))
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