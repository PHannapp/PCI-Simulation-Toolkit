import Equations
import Calculate_PCI  # Module for calculating Pressure-Composition Isotherms
import numpy as np  # For numerical operations
import matplotlib.pyplot as plt  # For plotting
import pandas as pd
import time

all_data = pd.DataFrame()

plt.rcParams.update({'font.size': 14})

# Initialize structure with given solids and their properties --------------------------------
solids = [
    ["CE", "LA"],
    ["AL", "NI"]
    ] # Define solids involved in the structure

molM = [
    [140.115, 138.9055], 
    [26.98, 58.6934]
    ] # g / mol

color = ["blue", "orange", "red", "green", "yellow", "purple", "pink", "black", "white", "gray"]
max_total = [1, 1.5/5]
max_concentrations = [
    [1], 
    [1.25/5/2, 1.2/5/2]
    ]

def transform_to_actual_concentrations(vars, max_concentrations, max_total):
    total = []
    delete = 0
    for tot in max_total:
        if tot == 1:
            total.append(1)
        else:
            total.append(vars[-1] * tot)
            delete += 1
    vars = vars[:-delete]
    actual_concentrations = []
    start_idx = 0
    max_concentrations_conc = np.concatenate(max_concentrations)
    for n, sublattice in enumerate(max_concentrations):
        # fill sublattice with proportional concentration
        sub_x = vars[start_idx : start_idx + len(sublattice)]
        # make sum(sub_x) = 1
        if n > 0:
            sub_x = [x / max(sum(sub_x),1e-3) for x in sub_x]
        sub_max = max_concentrations_conc[start_idx : start_idx + len(sublattice)]
        sublattice_concentration = [sub_x[i] * total[n] for i in range(len(sub_x))]
        for i, conc in enumerate(sublattice_concentration):
            if conc > sub_max[i]:
                sublattice_concentration[i] = sub_max[i]
        # if the sum of the sublattice is bigger than 1:
        if sum(sublattice_concentration) > 1:
            sublattice_concentration = [x / sum(sublattice_concentration) for x in sublattice_concentration]
        # append the last element so the sum equals 1
        rest = 1 - sum(sublattice_concentration)
        sublattice_concentration.append(rest)
        actual_concentrations.append(sublattice_concentration)
        start_idx += len(sublattice)
    return actual_concentrations
#for i in range(10):
#    bounds = [np.random.rand() for _ in sum(max_concentrations, [])]
#    bounds.append(np.random.rand())
#    #bounds = [0,0,0,0,0,0]
#    print(bounds)
#    site_fractions = transform_to_actual_concentrations(bounds, max_concentrations, max_total)
#    print(site_fractions)
for i in range(1):
    x_min = -1 # -1 for deactivation
    x_max = 1.45
    y_min = -1 # -1 for deactivation
    y_max = 300e5
    y_scale = 'log' # log and lin
    x_unit = 'WT' # HM, HF, WT
    site_fractions =  [
    [
      0,1
    ],
    [0,1]
  ]
 # Initial site fractions for hydrogen
    multiplicity = [1, 5, 1, 6]  # Define the multiplicity of sites in the structure

#MHK: Avbsorption: MH1: p1 1 atm, T1 ~ 20 °C
#: Desorption: MH1: p2: so hoch wie möglich, T2 ~ 100°C
# MH2: p1 bisschen kleiner als p2_MH1 at 20°C
# MH2_des: p2 so hoch wie möglich at 100°C
# capacity at MH1_p2-dp gleich capacoty of MH2 p1

    n = 0
    Initial_Structure = Equations.Structure(
        solids, multiplicity, site_fractions
    )  # Create a structure instance
    T = 313  # Set temperature (in Kelvin)
    # Perform Gibbs energy minimization to find equilibrium compositions --------------------------
    t0 = time.time()
    Calculate_PCI.gm, Calculate_PCI.x_H, fractions = Equations.gibbs_minimizer(
        T, Initial_Structure,
        #[40296.817431674805, -28564.60625096561, -118071.08470909395, -130056.37204219477, 0.0, -332349.51304057473, -78974.79537804768, 0.0, -270465.1659823443, 0.0, 140687.6162067526, 0.0, -20652.471029125947, 0.0, 28706.095021482757, 0.0, 0.0, 0.0, 11742.091809772604, 0.0],
        steps=100
    )  # Minimize Gibbs energy for the structure at temperature T
    print("Time: ", time.time()-t0)
    HV_fractions = [[i[0] for i in fractions], [i[1] for i in fractions]]
    # Plot Gibbs Energy vs. Hydrogen Mole Fraction -------------------------------------------------
    plt.figure(1)
    plt.plot(
        Calculate_PCI.x_H,
        Calculate_PCI.gm,
        label=f"{site_fractions}",
        color=color[n],
    )
    plt.ylabel("Gibbs Energy, $\t{\Delta G}$ / J/mol")
    plt.xlabel("Hydrogen mole fraction, $\t{x_{H}}$ / ")
    #plt.legend()
    plt.tight_layout()  # Adjust layout
    #plt.show()
    # Plot site_fractions vs. Hydrogen Mole Fraction -------------------------------------------------
    plt.figure(2)
    for i, site_fraction in enumerate(HV_fractions):
        if i == 0:
            plt.plot(
                Calculate_PCI.x_H,
                site_fraction,
                "--",
                label=f"{site_fractions}: 2b site",
                color=color[n],
            )
        else:
            plt.plot(
                Calculate_PCI.x_H,
                site_fraction,
                label=f"{site_fractions}: 6c site",
                color=color[n],
            )

    plt.ylabel("site occupancy of interstitial sites, $\t{y}$ / ")
    plt.xlabel("Hydrogen mole fraction, $\t{x_{H}}$ / ")
    plt.tight_layout()  # Adjust layout
    #plt.legend()
    # Calculate Plateau Pressures and related properties for Pressure-Composition Isotherms (PCIs) -
    BOUND_MAX = Calculate_PCI.x_H[-1]  # Get the maximum hydrogen mole fraction
    CTs, plateau_slopes, plateau_intercepts, plateau_pressures = Calculate_PCI.plateau(
        T
    )  # Calculate plateau pressures and related properties at temperature T
    print("Plateau pressure(s): ", plateau_pressures)
    print("CTs(s): ", CTs)
    # Plotting Pressure vs. Hydrogen Mole Fraction in a semilog plot -------------------------------
    if plateau_pressures:
        pl_root = np.log10(
            plateau_pressures[0]
        )  # Calculate the logarithmic base for plateau pressures
        p_plot = np.logspace(
            pl_root - 2, pl_root + 2, 50
        )  # Generate pressures for plotting
        # Adjust pressure points for detailed plotting around plateau pressures
        for k in plateau_pressures:
            p_plot = np.append(p_plot, k + k / 200000000)
            p_plot = np.append(p_plot, k - k / 200000000)
        p_plot = np.sort(p_plot)  # Sort the pressures for plotting
        x_plot = []  # Initialize list for hydrogen mole fractions
        # Calculate hydrogen mole fraction for each pressure point
        for p in p_plot:
            result, _, _ = Calculate_PCI.calculate_x_from_p(p, plateau_pressures, CTs, T)
            x_plot.append(Calculate_PCI.x_H[result])
        plt.figure(3)
        if x_unit == 'HM':
            x_plot = [- 6 * x / (x - 1) / 6 for x in x_plot]
            plt.xlabel("Hydrogen content, H$\t{per metal atom.}$ / ")
        if x_unit == 'HF':
            x_plot = [- 6 * x / (x - 1) for x in x_plot]
            plt.xlabel("Hydrogen content, H$\t{per mole f.u.}$ / ")
        if x_unit == 'WT':
            M_H = 1.008
            M_MH = 0
            for element in range(len(site_fractions[0])):
                M_MH += site_fractions[0][element] * molM[0][element]
            for element in range(len(site_fractions[1])):
                M_MH += site_fractions[1][element] * molM[1][element] * 5
            M_MH /= 6
            x_plot = [M_H * x / (M_H * x + M_MH * (1 - x)) * 100 for x in x_plot]
            plt.xlabel("Hydrogen content, $\t{wt.\%-}$H / ")
        if y_scale == 'log':
            plt.semilogy(
                x_plot, p_plot, label=f"{site_fractions}", color=color[n]
            )  # Create a semilogarithmic plot
        elif y_scale == 'lin':
            plt.plot(
                x_plot, p_plot, label=f"{site_fractions}", color=color[n]
            )  # Create a semilogarithmic plot
        if y_min != -1:
            plt.ylim(y_min, y_max)
        if x_min != -1:
            plt.xlim(x_min, x_max)
        plt.ylabel("Hydrogen Pressure, $\t{P}$ / Pa")
        plt.tight_layout()  # Adjust
            # Create a DataFrame for the current solid and site
        temp_df = pd.DataFrame({
            'x_H': x_plot,
            'p': p_plot
        })
        # Append to the main DataFrame
        all_data = pd.concat([all_data, temp_df], ignore_index=True)

        #plt.legend()
        #plt.grid()
        # Save the DataFrame to a CSV file
        all_data.to_csv('/home/pet69516/Dokumente/DFTforCALPHAD/CALPHAD_CALC/solid_site_data.csv', index=False)
    #plt.show()
    plt.figure(1)
    plt.savefig('/home/pet69516/Dokumente/DFTforCALPHAD/CALPHAD_CALC/plot_gibbs.png')  
    plt.figure(2)
    plt.savefig('/home/pet69516/Dokumente/DFTforCALPHAD/CALPHAD_CALC/plot_sf.png')
    plt.figure(3)
    plt.savefig('/home/pet69516/Dokumente/DFTforCALPHAD/CALPHAD_CALC/plot_PCI.png')
