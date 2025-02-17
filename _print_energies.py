from pymongo import MongoClient
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution  # Import optimization functions

# Connect to MongoDB
client = MongoClient('mongodb://localhost:27017/')
db = client.DFTforCALPHAD
collection_system = db.System
collection_calculation = db.calculation


def objective_function2(fit, x_values, y_values):
    RMSE = 0
    for x, y in zip(x_values, y_values):
        fit_y = x*(1-x)*(fit[0]+fit[1]*(x-(1-x))+fit[2]*(x-(1-x))**2+fit[3]*(x-(1-x))**3+fit[4]*(x-(1-x))**4+fit[5]*(x-(1-x))**5+fit[6]*(x-(1-x))**6+fit[7]*(x-(1-x))**7+fit[8]*(x-(1-x))**8+fit[9]*(x-(1-x))**9+fit[10]*(x-(1-x))**10+fit[11]*(x-(1-x))**11+fit[12]*(x-(1-x))**12+fit[13]*(x-(1-x))**13+fit[14]*(x-(1-x))**14)
        if fit_y > y:
            RMSE += (y - fit_y)**2
        else:
            RMSE += (y - fit_y)**2
    return (RMSE**(1/2)) / len(y_values)

def print_enegies(keyword1, name = None, n_atoms=1):
    energies_dict = {"formation_enthalpy": {}}
    energies_dict_full = {"formation_enthalpy": {}}
    tasks = collection_calculation.find(
        {
            'type': "manual",
            keyword1: {'$exists': True}
        }
    )
    CALCULATED = False
    for task in tasks:
        if task.get('name'):
            if name not in task['name'] or not task["VASP"].get("formation_enthalpy") or task['name'] == "POSCAR_FCC_0_1":
                continue
        print(task['name'], task["VASP"]["formation_enthalpy"])
        n_H_atoms = task['name'].split('_')
        if len(n_H_atoms) != 4:
            continue
        else:
            try:
                n_H_atoms = int(n_H_atoms[-2])
                #if n_H_atoms > 16:
                #    continue
            except ValueError as e:
                print(e)
                continue
            if not CALCULATED:
                base_atoms = task["number_of_atoms"] - n_H_atoms
                
        if energies_dict_full["formation_enthalpy"].get(n_H_atoms):
            energies_dict_full["formation_enthalpy"][n_H_atoms].append(task["VASP"]["formation_enthalpy_for_calphad"] * 1e3)
        else:
            energies_dict_full["formation_enthalpy"][n_H_atoms] = [task["VASP"]["formation_enthalpy_for_calphad"] * 1e3]
        if energies_dict["formation_enthalpy"].get(n_H_atoms):
            if energies_dict["formation_enthalpy"][n_H_atoms] > task["VASP"]["formation_enthalpy_for_calphad"] * 1e3:
                energies_dict["formation_enthalpy"][n_H_atoms] = task["VASP"]["formation_enthalpy_for_calphad"] * 1e3
        else:
            energies_dict["formation_enthalpy"][n_H_atoms] = task["VASP"]["formation_enthalpy_for_calphad"] * 1e3
    first_H = 6
    last_H = 29
    #for i in range(first_H):
    #    energies_dict["formation_enthalpy"][i] = energies_dict["formation_enthalpy"][first_H]
    #for i in range(last_H, 31+1):
    #    energies_dict["formation_enthalpy"][i] = energies_dict["formation_enthalpy"][last_H]
    energies_dict["formation_enthalpy"][0] = energies_dict["formation_enthalpy"][1]
    energies_dict_full["formation_enthalpy"][0] = [energies_dict["formation_enthalpy"][0]]
    # sort xs and also in the same order energies
    # Convert to sorted lists
    n_Hs, xs, energies = zip(*sorted(
        ((n_H, n_H / (n_H + base_atoms), energy) for n_H, energy in energies_dict["formation_enthalpy"].items()),
        key=lambda x: x[0]  # Sort by x_H values
    ))
    print(f'Energy if n_H = 0: {energies_dict["formation_enthalpy"][0]}')
    print(f'Energy if n_H = {n_Hs[-1]}: {energies_dict["formation_enthalpy"][n_Hs[-1]]}')
    energies_trimmed = []
    for ith_e, energy in enumerate(energies):
        n_H = n_Hs[ith_e]
        energy_trimmed = energy - n_H / n_Hs[-1] * energies_dict["formation_enthalpy"][n_Hs[-1]] - (1 - n_H / n_Hs[-1]) * energies_dict["formation_enthalpy"][0] 
        energies_trimmed.append(energy_trimmed)
    energies_full_trimmed = []
    for n_H in n_Hs:
        energy_trimmed_array = []
        for ith_e, energy in enumerate(energies_dict_full["formation_enthalpy"][n_H]):
            energy_trimmed = energy - n_H / n_Hs[-1] * energies_dict["formation_enthalpy"][n_Hs[-1]] - (1 - n_H / n_Hs[-1]) * energies_dict["formation_enthalpy"][0] 
            energy_trimmed_array.append(energy_trimmed)
        energies_full_trimmed.append(energy_trimmed_array)
    
    for n_H, energies_full in energies_dict_full["formation_enthalpy"].items():
        plt.plot([n_H / 30 for _ in energies_full], energies_full, 'o')    
    plt.xlabel("Hydrogen occupancy")
    plt.ylabel("J / mol (POSCAR)")
    plt.close()

    x_values = [n_H / max(n_Hs) for n_H in n_Hs]
    for ith_nH, n_H in enumerate(n_Hs):
        plt.plot([n_H / 31 for _ in energies_full_trimmed[ith_nH]], energies_full_trimmed[ith_nH], 'o', color='grey')    
    plt.plot(x_values, energies_trimmed, color='green', label="Minimum energy trajectory")
    #plt.plot(x_values, energies_trimmed, 'o')

    # Fitting equation: =x*(1-x)*(fit[0]+fit[1]*(x-(1-x))+fit[2]*(x-(1-x))**2)
    bounds = []
    bound = 1e7
    parameters = 3
    total = 15
    for i in range(parameters):
        bounds.append((-1*bound, bound))
    for i in range(total-parameters):
        bounds.append((0,0))
    result = differential_evolution(
        objective_function2, args=(x_values, energies_trimmed), bounds=bounds, workers=1, disp=True, popsize=8, maxiter=4000
    )
    print(result.fun)
    print(result.x)
    y_fit_values = [x*(1-x)*(result.x[0]+result.x[1]*(x-(1-x))+result.x[2]*(x-(1-x))**2) for x in np.linspace(0,1,100)]
    plt.plot(np.linspace(0,1,100), y_fit_values, color='red', label="Polynomial fit")
    plt.xlabel("Hydrogen occupancy")
    plt.ylabel("J / mol (POSCAR)")
    plt.legend()
    plt.savefig(r"/home/pet69516/Dokumente/DFTforCALPHAD_home/PLOTS/" + name + ".png")
    plt.show()

print_enegies('VASP', name='C14_', n_atoms=4)