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
"database_manual_C14_min",
#"database_TiMn2_no",
#"database_TiMn2_metal",
]
names = [
"DFT calculated",
#"No interactions",
#"Metal interaction",
]
multiplicitiess = [
 [[8, 16, 31]]  
]
plot_descriptor = "C14_manual_min"

exp_colors = ["black", "green"]
exp= [
[
[0.012363356180875473, 29.5645084808163],
[0.025773311716455027, 52.92206382714276],
[0.08523363885901411, 61.59610760786937],
[0.1167600365988731, 99.89582311955257],
[0.16007833443022937, 317.1935162531259],
[0.2625744417869238, 532.8517760256739],
[0.6802317928632197, 845.8939564350841],
[0.8387081239867088, 1001.449937500863],
[0.9443632919723258, 1134.2083555634283],
[0.9874440181710192, 1855.8077975584115],
[0.9268471997046407, 66.86348579352503],
[0.9085541839895982, 45.087272842991155],
[0.8787678379376213, 35.727493780216925],
[0.7054015442156123, 27.10532139640931],
[0.5353644637783521, 22.290154116834014],
[0.2844888196863412, 19.35674958666203],
[0.1341477117678219, 11.841579485610012],
[0.07780470969709616, 6.105076703933677]],[
[0.07437597316082636, 4.266408567687588],
[0.051549833860378536, 9.221076195004322],
[0.09813313642711527, 26.53043198888247],
[0.10505160762155483, 64.40675359581314],
[0.10680450101931069, 85.79030410575871],
[0.16145560781418042, 148.1083837512009],
[0.2887651090742732, 294.9138507924815],
[0.48854037915148396, 395.6603250891027],
[0.6882674928166685, 464.07082441414633],
[0.83850265662873, 564.420534802508],
[0.9227346421175978, 693.0574214368733],
[0.9839735460776601, 1154.305640974251],
[1.0220748992728383, 1742.4522700815946],
[0.8657078189960993, 52.99801395374173],
[0.7847536799524857, 40.53620375859101],
[0.5817326677047048, 35.18651298500826],
[0.3786988137470504, 29.467631707668694],
[0.21194921103744965, 23.379044052107623],
[0.08796571263463729, 12.613240745548286],
[0.059656163218132496, 6.161022463588862],
[0.0628473281217394, 4.543020070110604],
]
]
for ex in exp:
    for e in ex:
        print(e)
        e[1] = e[1]*1e5
x_Tis = [1]
Ts = [293]
exp_names = ["1st cylce - 293 K", "2nd cycle - 293 K"]#, "LaNi4.5Al0.5", "LaNi4Al"]
for num in range(len(file_names)):
    for xps in range(len(x_Tis)):
        file_name = file_names[num]

        database_modules = [file_name] 
        name = [names[num]]
        solids = [[["TI"], ["MN"]]]  # Define solids involved in the structure
        multiplicities = multiplicitiess[num]  # Define the multiplicity of sites in the structure
        x_Ti = x_Tis[xps]
        site_fractions = [[[1], [1]]]  # Initial site fractions for hydrogen
        T = Ts[xps]  # Set temperature (in Kelvin)
        color = ["blue"]
        FORCE = True
        MIN_Y = 4
        MAX_Y = 10
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