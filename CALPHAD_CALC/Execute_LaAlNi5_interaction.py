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
#"database_interaction_trimmed",
"database_interaction_metal",
#"database_interaction_metal_12",
#"database_interaction_metal_12_13",
]
names = [
"No interaction",
#"(H,Va) interaction trimmed",
"Metal interaction",
#"Metal interaction_12",
#"Metal interaction_12_13",
]
multiplicitiess = [
[[1,5,1,6]] ,
[[1,5,1,6]] ,
[[1,5,1,6]] 
]
plot_descriptor = "LaAlNi5_comb_313_both_exp"

exp_colors = ["red", "black", "green"]
exp=[]
exp2 = [
[
[0.0014418807752160923,  0.01897750129952208],
[0.022305192236116822,  0.1135568882786295],
[0.0380252094683505,  0.3173811188041542],
[0.09314784260054337,  0.7430718313581607],
[0.13457260836818374,  0.9828135943840031],
[0.2823917636956695,  1.2372547772894387],
[0.3494099132233948,  1.2541169201774856],
[0.5557600499383097,  1.445709394703004],
[0.6397730503105026,  1.6430103660956341],
[0.7790728003258887,  2.148375913159975],
[0.867297147596426,  2.9535029295251096],
[0.9075262073553505,  6.91361289798559],
[0.9222175657418505,  16.594007901977037],
[0.9181937805703635,  6.090702893079931],
[0.9119280466325344,  3.6664622519614993],
[0.8875919125727902,  2.0709959538221625],
[0.8174643412920073,  1.4322316221939881],
[0.737744746073659,  1.0686912873592294],
[0.509083555818663,  0.8163754725523986],
[0.41654528883040615,  0.7653003044150752],
[0.2963270138707759,  0.783770314039761],
[0.12826584530260454,  0.7066196630520573],
[0.09638328238883309,  0.6143332141396091],
[0.07835684204665022,  0.4765536931406123],
[0.05821740162533964,  0.3469443597626429],
[0.041339776860158134,  0.1862920793148639],
[0.01613323916171633,  0.04554967296705381],
[0.005726493973114227,  0.016717303142359544],
]] 
exp=[[
[0.0006330208281436767,  0.006293550629977665],
[0.012118245946175665,  0.01609431541844492],
[0.036328361970687645,  0.049164811947482304],
[0.05536587724671113,  0.0796389753040311],
[0.09780464860017413,  0.13068794197690756],
[0.18173559071686685,  0.2118674987874928],
[0.288065505933105,  0.26657764868803463],
[0.41673871186122785,  0.33128104552686394],
[0.5198947309808013,  0.3862608068429524],
[0.6825517772938945,  0.6189186573372221],
[0.7579955512702912,  0.9064388248818516],
[0.8227483068158172,  1.667854163808349],
[0.8704417371732689,  3.6183369027815955],
[0.9117580688175699,  7.652538270878649],
[0.9254383522702296,  14.617483733155996],
[0.9139472658482333,  5.862935178215428],
[0.8907101262817938,  2.8440516385051375],
[0.8546924134212137,  1.4697303259431247],
[0.8027348844297391,  0.7037074670125388],
[0.7390313022938213,  0.4074990258387027],
[0.5965840320496101,  0.2576289324879985],
[0.5444741091550638,  0.23858859281120384],
[0.5008689383127065,  0.2266578993064452],
[0.33073286814117536,  0.16262330456225718],
[0.2392906649942413,  0.13259425211811487],
[0.10432241860846786,  0.07291685074041433],
[0.05235902831302883,  0.03580978048226369],
[0.021701477927794616,  0.015495171344145046],
[0.0059668074356502565,  0.00590713516329867],
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

exp_names = ["LaNi4.5Al0.5"]#, "LaNi4.5Al0.5", "LaNi4Al"]
for num in range(len(file_names)):
    file_name = file_names[num]

    database_modules = [file_name] 
    name = [names[num]]
    solids =[[["CE", "LA"], ["AL", "NI"]]]  # Define solids involved in the structure
    multiplicities = multiplicitiess[num]  # Define the multiplicity of sites in the structure
    x_Ce = 0
    x_Al = 0.1
    site_fractions = [[[x_Ce, 1 - x_Ce], [x_Al, 1 - x_Al]]]  # Initial site fractions for hydrogen
    color = ["blue"]
    FORCE = True
    MIN_Y = 3
    MAX_Y = 9
    x_unit = "HM"

    vars = [[]]

    dir = '/home/pet69516/Dokumente/DFTforCALPHAD/CALPHAD_CALC/'
    T = 313  # Set temperature (in Kelvin)
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
        plt.plot(global_structure.x_H, global_structure.gm, label=f"{global_structure.name[0]}", color=cmap(num))

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
plt.legend(loc='lower right', ncol=2)
plt.grid()
plt.savefig(f'{dir}/plot_PCI_{plot_descriptor}.png')
#plt.show()
plt.close()