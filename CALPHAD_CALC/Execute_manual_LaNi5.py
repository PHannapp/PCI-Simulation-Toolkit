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
"database_manual_LaNi5_min",
#"database_TiMn2_no",
#"database_TiMn2_metal",
]
names = [
"DFT calculated",
#"No interactions",
#"Metal interaction",
]
multiplicitiess = [
 [[4, 20, 31]]
]
plot_descriptor = "LaNi5_manual_min"

exp_colors = ["red", "blue"]
exp= [
[
[0.027136825288143346, 1617.2014493342422],
[0.02701890524553996, 2348.4767251400476],
[0.03845714937806688, 14277.061237921143],
[0.05565445623644871, 33538.15742273175],
[0.08813952603750619, 77813.20297193078],
[0.14274411350755067, 73172.87321194167],
[0.2060633725132185, 73061.44117121259],
[0.6056107117045153, 76850.0555538021],
[0.9025752215755642, 71846.83196074891],
[0.9527749248735216, 76210.58434017196],
[1.0214462322644455, 106570.97435072639],
[1.0463463806154665, 206459.6097652833],
[1.0592300962379704, 409835.7476969317],
[1.0280003043097874, 105279.53210568204],
[0.9801551979915557, 57745.7501416845],
[0.8677431625394654, 51958.515996844144],
[0.7836699760356044, 53977.77701153787],
[0.3699227813914566, 52584.82536073639],
[0.15266080870326, 54803.829876196054],
[0.07849290577808207, 44208.88164821543],
[0.06218190117539654, 36043.74037563902],
[0.052504849937236114, 22547.492649763582],
[0.03321921716307202, 7104.895354235898]],[
[0.03258777435429269, 1656.3804844640251],
[0.04331849823119935, 94437.72102421409],
[0.05724067100308114, 221860.56239643582],
[0.08325915782266347, 395074.70040968957],
[0.12683251550077984, 533194.6335478964],
[0.1497318269998859, 579738.0026732252],
[0.32003879949788894, 577366.4681935669],
[0.48270759633306703, 568230.7161577771],
[0.687934116931036, 593315.1043354771],
[0.9509985164897867, 664967.486146311],
[1.0098406177488686, 941354.6686683025],
[1.0477766366160755, 2236986.4697075016],
[1.0250446955000192, 1211575.4111192285],
[0.9826771653543307, 625657.811096677],
[0.8386853817185895, 437570.04788439017],
[0.5723040054775762, 445716.32557396],
[0.3626992278139147, 442611.8646752856],
[0.16073262581307773, 444768.68346851954],
[0.11169690745178595, 333588.6217511441],
[0.06166457453687854, 185198.15428330185],
[0.049899197382935856, 85756.06872168115],
[0.04029441971927422, 42680.48771370625],
[0.04039332040016738, 31213.310884331575],

]
]

x_Tis = [1, 1]
Ts = [273, 323]
exp_names = ["LaNi5 - 273 K", "LaNi5 - 323 K"]#, "LaNi4.5Al0.5", "LaNi4Al"]
for num in range(len(file_names)):
    for xps in range(len(x_Tis)):
        file_name = file_names[num]

        database_modules = [file_name] 
        name = [names[num]]
        solids = [[["LA"], ["NI"]]]  # Define solids involved in the structure
        multiplicities = multiplicitiess[num]  # Define the multiplicity of sites in the structure
        x_Ti = x_Tis[xps]
        site_fractions = [[[1], [1]]]  # Initial site fractions for hydrogen
        T = Ts[xps]  # Set temperature (in Kelvin)
        color = ["blue"]
        FORCE = True
        MIN_Y = 3
        MAX_Y = 8
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