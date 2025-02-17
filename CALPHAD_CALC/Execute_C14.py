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
"database_TiMn2_first",
#"database_TiMn2_no",
#"database_TiMn2_metal",
]
names = [
"Normal",
#"No interactions",
#"Metal interaction",
]
multiplicitiess = [
 [[1,2,1,1,2]]  ,
 [[1,2,1,1,2]]  ,
 [[1,2,1,1,2]]  ,
]
plot_descriptor = "Ti02Zr08Mn2_440_interaction"

exp_colors = ["red", "black", "green"]
exp = [
[
[0.13746782440457295, 13.331792202902768],
[0.33074298712670347, 120.56298229783027],
[0.9644757521911925, 441.7737677762894],
[1.9175838235581741, 633.2494682845097],
[2.6599485862537766, 1040.9309615315826],
[2.7580984387094403, 1897.3191822677054],
[2.4044282195173845, 150.02485917480755],
[1.779247174572038, 95.2143338768725],
[0.7871658963679619, 72.3913545618711],
[0.23609556710218882, 41.8489712577209],
[0.11405120254891965, 13.525440345301943],
]
,
[[0.04403925408978204, 0.683120297239543],
[0.04234431614039802, 1.0059549596322488],
[0.15374265171024626, 1.482022179685194],
[0.7063800924136554, 2.1013720168190115],
[1.2935299765346906, 2.6076574892314643],
[2.1962889960582332, 4.119741817857355],
[2.7484004215359716, 6.586905271808742],
[3.231428514042273, 13.934456110634986],
[3.33122529157927, 21.944528275473395],
[3.4529393701858844, 40.56508060881823],

]]
tmp = [[
[0.001008144281858414,  0.0012406910156891433],
[0.013615809109052496,  0.0024618056394532717],
[0.03893957288678018,  0.0060614686049115804],
[0.09729764580726274,  0.011732775738179532],
[0.2408292572848682,  0.016981310056119473],
[0.3907614126914815,  0.02277818035163065],
[0.5012411311144391,  0.04525310591898317],
[0.5848497015130956,  0.2961665747411514],
[0.6577438082650249,  2.695351124920116],
[0.7126026827188245,  19.766572780110028],
[0.6738741167747591,  1.2753675145827967],
[0.6159527109996161,  0.09951263084539025],
[0.5704426163688636,  0.03604591994030118],
[0.4652966845534126,  0.017029736446005715],
[0.34406440400795973,  0.01405708763337101],
[0.21643157953349884,  0.012520029848664958],
[0.11327262976194319,  0.010875041796954563],
[0.04961300740575755,  0.005206201397570423],
[0.017894561002986364,  0.0022243350311715528],
[0.005269312263899373,  0.0012096725650480233],
]
]
for expe in exp:
    for i in range(len(expe)):
        expe[i][1] *= 1e5
        expe[i][0] /= 3
x_Tis = [1,0]
Ts = [441, 317]
exp_names = ["TiMn2 - 441 K", "Ti0.2Zr0.8Mn2 - 317 K"]#, "LaNi4.5Al0.5", "LaNi4Al"]
for num in range(len(file_names)):
    for xps in range(len(x_Tis)):
        file_name = file_names[num]

        database_modules = [file_name] 
        name = [names[num]]
        solids = [[["TI", "ZR"], ["MN"]]]  # Define solids involved in the structure
        multiplicities = multiplicitiess[num]  # Define the multiplicity of sites in the structure
        x_Ti = x_Tis[xps]
        site_fractions = [[[x_Ti, 1 - x_Ti], [1]]]  # Initial site fractions for hydrogen
        T = Ts[xps]  # Set temperature (in Kelvin)
        color = ["blue"]
        FORCE = True
        MIN_Y = 3
        MAX_Y = 11
        x_unit = "HM"

        vars = [[]]

        dir = '/home/pet69516/Dokumente/DFTforCALPHAD/CALPHAD_CALC/'
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
    plt.ylim(np.power(10,MIN_Y), np.power(10,MAX_Y))
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