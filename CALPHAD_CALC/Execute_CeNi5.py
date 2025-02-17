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
"database_interaction_no",
#"database_interaction_no_Entropy",
#"database_interaction_trimmed",
#"database_interaction_trimmed_Entropy"
#"database_VASP_prec",
#"database_VASP_GGA=91_eff",
#"database_VASP_GGA=91_prec",
#"database_VASP_METAGGA=R2SCAN_eff",
#"database_VASP_METAGGA=R2SCAN_prec",
#"database_VASP_eff_3SL",
#"database_VASP_prec_3SL"
]
names = [
"PBE",
#"PBE_Entropy",
#"PBE_trimmed",
#"PBE_trimmed_Entropy",
#"PBE (HQ)",
#"PW91",
#"PW91 (HQ)",
#"R2Scan",
#"R2SCAN (HQ)",
#"PBE, 3SL Model",
#"PBE (HQ), 3SL Model"
]

multiplicitiess = [
[[1,5,1,6]] ,
[[1,5,1,6]] ,
[[1,5,1,6]] ,
[[1,5,1,6]] ,
[[1,5,1,6]] ,
[[1,5,1,6]] ,
#[[1,5,1,3,3]] ,
#[[1,5,1,3,3]] 
]
plot_descriptor = "CeNi5_comb_296_no_exp"

exp_colors = ["black", "red"]
exp=[]
exp2 = [
[
[0.0035113547065084294, 986779.0572165092],
[0.030435885448269695, 14547382.028839834],
[0.0821960390327362, 28617685.21719321],
[0.17018631069023488, 45403578.488040656],
[0.24517561747122774, 64861588.03311403],
[0.41652326148153307, 68354422.86056784],
[0.5889924500900219, 69803667.00512294],
[0.7177835692472968, 72035348.90659674],
[0.9077200067640826, 73562635.70436288],
[1.0048616844555411, 75122303.89540397],
[1.0516308402383343, 90260045.92237107],
[1.0635301548775997, 101299999.9999999],
[1.0607673255015864, 59018110.28048984],
[1.0409848703385025, 21446170.98887003],
[1.024522286657847, 7433837.955410418],
[0.9805084998657129, 3384756.577031162],
[0.7536282340770509, 2758602.7957290984],
[0.5299236056539773, 2496945.077566524],
[0.3388705971292437, 2510076.1321984404],
[0.2832733186778208, 2308024.554498622],
[0.22121534651003172, 1931041.8130569041],
[0.20634431170484724, 1249469.3825423268],
[0.20171638599039105, 520374.00765385997]
],
[
[0.19627776506749162, 509570.14854548837],
[0.2062771682366633, 4241088.420006282],
[0.21820881121246183, 14547382.028839834],
[0.29860689737493906, 21900870.72047555],
[0.4404835324427291, 22959504.99800376],
[0.6402055087485452, 24451042.046124797],
[0.7613397857377325, 25632945.33945418],
[0.9751395092061154, 30798187.801669735],
[1.0086689677810825, 42858121.88078572],
[1.0635450756483071, 99718498.92139088],
[1.0275163879798273, 31616565.040502466],
[0.9924500900219834, 11488918.51529677],
[0.9053973401239418, 8519980.534409262],
[0.7079085058340212, 7591449.499442121],
[0.6478747848922223, 7512230.389540397],
[0.4721404342938994, 7279498.6980235465],
[0.2943122022062846, 6418473.741701322],
[0.22112084829555062, 2133397.8836891903],
[0.20735146372760643, 1365997.77534616],
[0.1933184788771623, 1154933.6991474503],
[0.20391968646487146, 509570.14854548837]
]
]

exp_names = ["Experiment, 1st cycle", "Experiment, activated"]
for num in range(len(file_names)):
    file_name = file_names[num]

    database_modules = [file_name] 
    name = [names[num]]
    solids =[[["CE", "LA"], ["AL", "NI"]]]  # Define solids involved in the structure
    multiplicities = multiplicitiess[num]  # Define the multiplicity of sites in the structure
    x_Ce = 1
    x_Al = 0
    site_fractions = [[[x_Ce, 1 - x_Ce], [x_Al, 1 - x_Al]]]  # Initial site fractions for hydrogen
    color = ["blue"]
    FORCE = True
    MIN_Y = 5
    MAX_Y = 10
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

        plt.figure(2)
        plt.plot(global_structure.x_H, global_structure.gm,
            label=f"{global_structure.name[0]}", color=cmap(num))

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
plt.legend(loc='upper left', ncol=2)
plt.grid()
plt.savefig(f'{dir}/plot_PCI_{plot_descriptor}.png')
#plt.show()
plt.close()