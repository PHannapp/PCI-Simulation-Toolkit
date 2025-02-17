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
"database_manual_BCC2",
#"database_TiMn2_no",
#"database_TiMn2_metal",
]
names = [
"DFT calculated",
#"No interactions",
#"Metal interaction",
]
multiplicitiess = [
 [[16, 31]]  
]
plot_descriptor = "BCC_manual"

exp_colors = ["red", "blue", "green"]
exp= [
[
[0.000840108401084014,0.00015082254718489071],
[0.010840108401083959,0.0005365890883515254],
[0.005420054200541979,0.0008719600631335286],
[0.03523035230352306,0.0013641102857447963],
[0.09756097560975607,0.001362523809748111],
[0.19241192411924118,0.0013601131466291238],
[0.2818428184281842,0.001169459155437193],
[0.3631436314363144,0.0013060942307295668],
[0.43902439024390255,0.001631761037427775],
[0.46070460704607047,0.004640179996277237],
[0.8346883468834692,3172.052486339696],
[0.9132791327913281,48360.07057131725],
[0.9837398373983743,207180.10811823385],
[1.4769647696476969,229613.5869010556],
[1.7425474254742552,212037.06182872143],
[1.9349593495934962,203533.29560211924],
[1.9674796747967485,219181.66356349696]],[
[0.00,0.05711586478126423],
[0.02439024390243899,0.2272788367979062],
[0.027100271002710008,0.6966642282662012],
[0.051490514905149,1.8384311677807883],
[0.11653116531165308,3.4641420575026345],
[0.1788617886178862,4.170339301617164],
[0.3821138211382114,3.855591099246404],
[0.4254742547425476,5.596307377148788],
[0.4769647696476966,28.906637944294758],
[0.4959349593495935,52.51740156005165],
[0.7479674796747968,38773.052517637785],
[0.9403794037940383,7756621.495413083],
[1.0054200542005423,12587980.272093719],
[1.384823848238483,12974661.100517929],
[1.9701897018970198,11910154.103838412]],[
[0.005420054200542035,227.4053632687331],
[0.008130081300812997,4508.0438777575155],
[0.05420054200542007,43934.306895163325],
[0.1436314363143632,181257.03863688442],
[0.2628726287262873,597374.910808231],
[0.5067750677506775,8425986.184204623],
[0.6910569105691058,95100022.49027129],

]
]

x_Tis = [1,1,1]
Ts = [298, 393, 823]
exp_names = ["V - 298 K", "V - 393 K", "V - 823 K"]#, "LaNi4.5Al0.5", "LaNi4Al"]
for num in range(len(file_names)):
    for xps in range(len(x_Tis)):
        file_name = file_names[num]

        database_modules = [file_name] 
        name = [names[num]]
        solids = [[["VV"]]]  # Define solids involved in the structure
        multiplicities = multiplicitiess[num]  # Define the multiplicity of sites in the structure
        x_Ti = x_Tis[xps]
        site_fractions = [[[1]]]  # Initial site fractions for hydrogen
        T = Ts[xps]  # Set temperature (in Kelvin)
        color = ["blue"]
        FORCE = True
        MIN_Y = -6
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