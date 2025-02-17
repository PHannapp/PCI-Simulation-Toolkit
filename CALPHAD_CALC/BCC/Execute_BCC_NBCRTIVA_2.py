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
#"database_BCC_1_nothing",
"database_BCC_1",
#"database_BCC_1_both",
]
file_names_FCC = [
#"database_FCC_2_nothing",
"database_FCC_2",
#"database_FCC_2_both",
]
names_BCC = [
"Simulation",
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
exp= [
[
[0.8436131611880859, 0.24067529130834098],
[0.9163664794267892, 3.584131687611682],
[0.9355609718982769, 10.400045593802027],
[0.9585177528066309, 21.423730898274894],
[1.0347657490076116, 32.51495352988446],
[1.1642944659728038, 33.693147097384966],
[1.3319340628779672, 39.10856036062965],
[1.4996103906325846, 59.25216947079542],
[1.5529855621525335, 80.2591478219607],
[1.629181085711438, 83.2518234784091],
[1.6863237929324604, 83.16134427268078],
[1.6291286130693616, 56.899039098749405],
[1.6138170961115157, 33.40615677696254],
[1.590865562467369, 16.845951107433788],
[1.5450884295200593, 10.678905709117306],
[1.499353274686411, 9.178865918124911],
[1.4269305340927874, 6.778852058162478],
[1.3240421775097013, 5.405463252125752],
[1.1868786911224163, 4.654261701350002],
[1.0459109381846043, 4.163216970347567],
[0.9658586754330964, 2.849714668914963],
[0.923901550828937, 1.9492110953972157],
[0.8856332529627369, 0.5555467308434392],
[0.8740525408565114, 0.18427968520027602],
[0.809091409966129, 0.04344253473005345],
[0.7823251152430404, 0.02109049530981154],
[0.6906134314221926, 0.0027057193254737735],



]
]
for expe in exp:
    for i in range(len(expe)):
        expe[i][1] *= 1.01325e5

exp_names = ["Nb55.4Cr15Ti28.3V1.3, Experiment"]
for num in range(len(file_names_BCC)):
    file_name_BCC = file_names_BCC[num]
    file_name_FCC = file_names_FCC[num]

    database_modules = [file_name_BCC, file_name_FCC] 
    name = [names_BCC[num], names_FCC[num]]
    solids =[[["NB", "CR", "TI", "VA"]], [["NB", "CR", "TI", "VA"]]]  # Define solids involved in the structure
    multiplicities = [multiplicitiess_BCC[num], multiplicitiess_FCC[num]]  # Define the multiplicity of sites in the structure
    x_NB = 0.554
    x_CR = 0.15
    x_TI = 0.283  
    site_fractions = [[[x_NB, x_CR, x_TI, 1-x_NB-x_CR-x_TI]], [[x_NB, x_CR, x_TI, 1-x_NB-x_CR-x_TI]]]  # Initial site fractions for hydrogen
    color = ["blue"]
    FORCE = True
    MIN_Y = -1
    MAX_Y = 8
    MIN_X = 0
    MAX_X = 2
    x_unit = "HM"

    vars = [[], []]

    T = 273.15+200  # Set temperature (in Kelvin)
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