from pymongo import MongoClient
from bson.objectid import ObjectId  # Required to work with _id
from write import clean_mcsqs_poscar
import math
from scipy.optimize import minimize, differential_evolution  # Import optimization functions
import matplotlib.pyplot as plt
import numpy as np
from CONSTANTS import AVOGADRO, EV_TO_KJ

# Connect to MongoDB
client = MongoClient('mongodb://localhost:27017/')
db = client.DFTforCALPHAD
collection_system = db.System
collection_calculation = db.calculation

def get_nested_value(dictionary, keys):
    """
    A helper function to traverse a dictionary with a list of keys.
    It will return the value if all keys are found, or None otherwise.
    """
    value = dictionary
    for key in keys:
        if isinstance(value, dict) and key in value:
            value = value[key]
        else:
            return None
    return value

# Function to compute GCD of a list of numbers
def find_gcd(numbers):
    gcd_result = numbers[0]
    for num in numbers[1:]:
        gcd_result = math.gcd(gcd_result, num)
    return gcd_result

def calculate_formation_enthalpy(task, source, FORCE_reference=None):
    poscar_content = task["POSCAR"]
    if not poscar_content:
        return None, None, None
    # Split the content into individual lines
    poscar_lines = poscar_content.splitlines()
    
    # Extract the elements (line 5) and the number of atoms (line 6) from the POSCAR
    elements = poscar_lines[5].split()
    number_of_atoms = [int(i) for i in poscar_lines[6].split()]
    total_number_of_atoms = sum(number_of_atoms)
    # dH = E / at
    # Split the source key on "." to handle nested dictionaries
    source_parts = source.split(".")
    # Try to get the formation enthalpy using the nested source keys
    formation_enthalpy = get_nested_value(task, source_parts)
    if formation_enthalpy is None or "energy" not in formation_enthalpy:
        raise KeyError(f"Could not find energy for the source {source}")
    if formation_enthalpy["energy"] is None:
        return None, None, None
    formation_enthalpy = formation_enthalpy["energy"] / total_number_of_atoms
    for site in range(len(number_of_atoms)):
        # check if reference element is available
        reference_element = collection_calculation.find_one({
            'type': 'reference',
            #'system_id': task['system_id'],
            'elements': elements[site]
        })
        if reference_element:
            if FORCE_reference:
                ref_source_parts = FORCE_reference.split(".")
            else:
                ref_source_parts = source.split(".")
            ref_source = get_nested_value(reference_element, ref_source_parts)
            if ref_source is None:
                # Fallback to "VASP" if specific source key isn't available in the reference
                ref_source = reference_element.get("VASP", None)
            
            if ref_source and "energy" in ref_source:
                reference_energy = ref_source['energy']
                if 'VASP' in source_parts:
                    number_of_atoms_in_reference_unit_cell = int(reference_element['POSCAR'].splitlines()[6].split()[0])
                elif source == 'CHGNet_results':
                    number_of_atoms_in_reference_unit_cell = 1
                else:
                    print('Unknown source')
                    return None, None, None
                formation_enthalpy -= number_of_atoms[site] / total_number_of_atoms * reference_energy / number_of_atoms_in_reference_unit_cell
        else:
            return None, None, None
    if task['type'] == 'endmember':
        number_of_atoms_per_formula_unit = total_number_of_atoms / find_gcd(number_of_atoms)
    elif task['type'] == 'mcsqs':
        reference_endmembers = [[],[]]
        for elements in task['elements']:
            if len(elements) == 2:
                for i in range(2):
                    reference_endmembers[i].append([elements[i]])
            else:
                for i in range(2):
                    reference_endmembers[i].append(elements)
        number_of_atoms_per_formula_unit = 0
        count = 0
        for reference_endmember in reference_endmembers:
            endmember = collection_calculation.find_one({
                'type': 'endmember',
                'system_id': task['system_id'],
                'elements': reference_endmember
            })
            if endmember:
                # Get the nested value using the source parts
                nested_source_value = get_nested_value(endmember, source_parts)

                # Check if the nested value exists and contains the required key
                if nested_source_value and "number_of_atoms_per_formula_unit" in nested_source_value:
                    number_of_atoms_per_formula_unit += nested_source_value["number_of_atoms_per_formula_unit"]
                    count += 1
                    continue
        if not number_of_atoms_per_formula_unit:
            print('Unknown number_of_atoms_per_formula_unit')
            return None, None, None
        number_of_atoms_per_formula_unit /= count
    else:
        number_of_atoms_per_formula_unit = total_number_of_atoms
        print('Unknown type: ', task['type'])
    formation_enthalpy_calphad = formation_enthalpy * number_of_atoms_per_formula_unit * AVOGADRO * EV_TO_KJ

    return formation_enthalpy, formation_enthalpy_calphad, number_of_atoms_per_formula_unit


def calculate_mixing_enthalpy(task, source, FORCE_mcsqs=None):
    reference_endmembers = [[], []]
    for elements in task['elements']:
        if len(elements) == 2:
            for i in range(2):
                reference_endmembers[i].append([elements[i]])  # Keep as [elements[i]]
        else:
            for i in range(2):
                reference_endmembers[i].append(elements)  # Keep as elements

    formation_enthalpies, formation_enthalpies_calphad = [], []
    
    # Split the source key to handle nested keys
    source_parts = source.split(".")
    for reference_endmember in reference_endmembers:
        # Find the endmember in the collection
        endmember = collection_calculation.find_one(
            {
                'type': 'endmember',
                'system_id': task['system_id'],
                'elements': reference_endmember
            },
            sort=[(f'{source}.energy', 1)]  # Sort by 'VASP.eff.energy' in ascending order
        )
        
        # Safely get the nested values
        nested_source_value = get_nested_value(endmember, source_parts) if endmember else None
        
        if nested_source_value and 'formation_enthalpy' in nested_source_value:
            formation_enthalpies.append(nested_source_value["formation_enthalpy"])
            formation_enthalpies_calphad.append(nested_source_value.get("formation_enthalpy_for_calphad", 0))
        else:
            print(f'Formation_enthalpy from endmember {reference_endmember} missing.')
            return None, None

    # Safely access formation_enthalpy in task using the source_parts
    if FORCE_mcsqs:
        source_parts = FORCE_mcsqs.split(".")
    task_source_value = get_nested_value(task, source_parts)
    if not task_source_value or "formation_enthalpy" not in task_source_value:
        print(f'Formation enthalpy not found in task for source {source}.')
        return None, None

    mixing_enthalpy = task_source_value["formation_enthalpy"]
    
    # Calculate mixing enthalpy
    for proportion, formation_enthalpy in zip([task['mixing_proportion'], 1 - task['mixing_proportion']], formation_enthalpies):
        mixing_enthalpy -= proportion * formation_enthalpy
    
    # Calculate mixing enthalpy for CALPHAD
    mixing_enthalpy_calphad = mixing_enthalpy * task_source_value['number_of_atoms_per_formula_unit'] * AVOGADRO * EV_TO_KJ

    return mixing_enthalpy, mixing_enthalpy_calphad

def calculate_mixing_enthalpy_HV(task, source, FORCE_mcsqs=None):

    formation_enthalpies, formation_enthalpies_calphad = [], []
    
    # Split the source key to handle nested keys
    source_parts = source.split(".")
    for mixing_proportion in [0.125, 0.875]:
        # Find the endmember in the collection
        endmember = collection_calculation.find_one({
            'type': 'mcsqs',
            'system_id': task['system_id'],
            'elements': task['elements'],
            'mixing_proportion': mixing_proportion
        })
        
        # Safely get the nested values
        nested_source_value = get_nested_value(endmember, source_parts) if endmember else None
        
        if nested_source_value and 'formation_enthalpy' in nested_source_value:
            formation_enthalpies.append(nested_source_value["formation_enthalpy"])
            formation_enthalpies_calphad.append(nested_source_value.get("formation_enthalpy_for_calphad", 0))
        else:
            print(f'Formation_enthalpy from mixing_proportion {mixing_proportion} missing.')
            return None, None

    # Safely access formation_enthalpy in task using the source_parts
    if FORCE_mcsqs:
        source_parts = FORCE_mcsqs.split(".")
    task_source_value = get_nested_value(task, source_parts)
    if not task_source_value or "formation_enthalpy" not in task_source_value:
        print(f'Formation enthalpy not found in task for source {source}.')
        return None, None

    mixing_enthalpy = task_source_value["formation_enthalpy"]
    bounds = [0.125, 0.875]
    max_range = bounds[1] - bounds[0]
    # Calculate mixing enthalpy
    for proportion, formation_enthalpy in zip([(task['mixing_proportion'] - bounds[0]) / max_range, (bounds[1] - task['mixing_proportion']) / max_range], reversed(formation_enthalpies)):
        mixing_enthalpy -= proportion * formation_enthalpy

    # Calculate mixing enthalpy for CALPHAD
    mixing_enthalpy_calphad = mixing_enthalpy * task_source_value['number_of_atoms_per_formula_unit'] * AVOGADRO * EV_TO_KJ

    return mixing_enthalpy, mixing_enthalpy_calphad

def objective_function(fit, x_values, y_values):
    RMSE = 0
    for x, y in zip(x_values, y_values):
        fit_y = x*(1-x)*(fit[0]+fit[1]*(x-(1-x))+fit[2]*(x-(1-x))**2) + fit[3]
        RMSE += (y - fit_y)**2
    return (RMSE**(1/2)) / len(y_values)

def fit_interaction_parameters(x_values, y_values, key, bound, parameter_space, save):
    print(key)
    # Fitting equation: =x*(1-x)*(fit[0]+fit[1]*(x-(1-x))+fit[2]*(x-(1-x))**2)
    bounds = []
    for i in range(3):
        if parameter_space[i] == 1:
            bounds.append((-1*bound, bound))
        else:
            bounds.append((0,0))
    bounds.append((-1*bound, bound))
    result = differential_evolution(
        objective_function, args=(x_values, y_values), bounds=bounds, workers=1, popsize=8, maxiter=4000
    )
    print(result.x)
    y_fit_values = [x*(1-x)*(result.x[0]+result.x[1]*(x-(1-x))+result.x[2]*(x-(1-x))**2) for x in np.linspace(0,1,100)]
    y_values = [y - result.x[3] for y in y_values]
    if save is None:
        plt.plot(np.linspace(0,1,100), y_fit_values)
        plt.plot(x_values, y_values, 'o')
        plt.show()
    return result.x, np.linspace(0,1,100), y_fit_values

def objective_function2(fit, x_values, y_values):
    RMSE = 0
    for x, y in zip(x_values, y_values):
        fit_y = x*(1-x)*(fit[0]+fit[1]*(x-(1-x))+fit[2]*(x-(1-x))**2)
        RMSE += (y - fit_y)**2
    return (RMSE**(1/2)) / len(y_values)

def fit_interaction_parameters2(x_values, y_values, key, bound, parameter_space, save):
    print(key)
    # Fitting equation: =x*(1-x)*(fit[0]+fit[1]*(x-(1-x))+fit[2]*(x-(1-x))**2)
    bounds = []
    for i in range(3):
        if parameter_space[i] == 1:
            bounds.append((-1*bound, bound))
        else:
            bounds.append((0,0))
    result = differential_evolution(
        objective_function2, args=(x_values, y_values), bounds=bounds, workers=1, popsize=8, maxiter=4000
    )
    print(result.x)
    y_fit_values = [x*(1-x)*(result.x[0]+result.x[1]*(x-(1-x))+result.x[2]*(x-(1-x))**2) for x in np.linspace(0,1,100)]
    if save is None:
        plt.plot(np.linspace(0,1,100), y_fit_values)
        plt.plot(x_values, y_values, 'o')
        plt.show()
    return result.x, np.linspace(0,1,100), y_fit_values

def calculate_energies(sources, forces_thr, FORCE_mcsqs=None):
    for source in sources:
        tasks = collection_calculation.find(
            {
                f'{source}.energy': {'$exists': True},
                '$or': [
                    {f'{source}.formation_enthalpy': {'$exists': False}}, 
                    {f'{source}.mixing_enthalpy': {'$exists': True}}
                ],
                #'system_id': ObjectId('672c67ffafd1d7bf02525ca6'),
                #'_id': ObjectId('672da8624d88bce903cd3002'),
                #f'{source}.formation_enthalpy_for_calphad': {'$exists': False},
                'type' : {'$ne':'reference'}
            }
        )
        for task in tasks:
            if str(task['_id']) == '672da8624d88bce903cd3002':
                pass
            print("Found task: ", task['_id'])
            # Split the source key on "." to handle nested dictionaries
            source_parts = source.split(".")
            # Try to get the formation enthalpy using the nested source keys
            task_source = get_nested_value(task, source_parts)
            if task_source.get("total_forces"):
                if task_source["total_forces"] > forces_thr:
                    print("TOTAL FORCES TOO HIGH: ", task_source["total_forces"], task["elements"])                
                    collection_calculation.update_one(
                        {'_id':ObjectId(task["_id"])}, 
                        {
                            '$unset': {f"{source}.formation_enthalpy": "",
                                     f"{source}.formation_enthalpy_for_calphad": "",
                                     f"{source}.number_of_atoms_per_formula_unit": ""},

                            '$set': {f"{source}.FORCES_ERROR": True}
                        })
                    continue
            formation_enthalpy, formation_enthalpy_for_calphad, number_of_atoms_per_formula_unit = calculate_formation_enthalpy(task, source, FORCE_mcsqs)
            if formation_enthalpy:
                collection_calculation.update_one(
                    {'_id':ObjectId(task["_id"])}, 
                    {
                        '$set': {f"{source}.formation_enthalpy": formation_enthalpy,
                                 f"{source}.formation_enthalpy_for_calphad": formation_enthalpy_for_calphad,
                                 f"{source}.number_of_atoms_per_formula_unit": number_of_atoms_per_formula_unit}
                    })
            else:
                print("No formation enthalpy")
                #formation_enthalpy, formation_enthalpy_for_calphad, number_of_atoms_per_formula_unit = calculate_formation_enthalpy(task, source)

        tasks = collection_calculation.find(
            {
                f'{source}.formation_enthalpy_for_calphad': {'$exists': True},
                #'system_id': ObjectId('672c67ffafd1d7bf02525ca6'),
                #'_id': ObjectId('672da8624d88bce903cd3002'),
                #f'{source}.mixing_enthalpy_for_calphad': {'$exists': False},
                'type': 'mcsqs'
            }
        )
        if FORCE_mcsqs:
            tasks = collection_calculation.find(
                {
                    f'{FORCE_mcsqs}.formation_enthalpy_for_calphad': {'$exists': True},
                    #f'{source}.mixing_enthalpy_for_calphad': {'$exists': False},
                    'type': 'mcsqs'
                }
            )
        for task in tasks:
            formation_enthalpy = None
            for site in task['elements']:
                if site == ["H", "Va"]:
                    formation_enthalpy, formation_enthalpy_for_calphad = calculate_mixing_enthalpy(task, source, FORCE_mcsqs)
                    break
            if not formation_enthalpy and formation_enthalpy != 0:
                formation_enthalpy, formation_enthalpy_for_calphad = calculate_mixing_enthalpy(task, source, FORCE_mcsqs)
            if formation_enthalpy or formation_enthalpy == 0:
                if FORCE_mcsqs:
                    collection_calculation.update_one(
                        {'_id':ObjectId(task["_id"])}, 
                        {
                            '$set': {f"{FORCE_mcsqs}.mixing_enthalpy": formation_enthalpy,
                                     f"{FORCE_mcsqs}.mixing_enthalpy_for_calphad": formation_enthalpy_for_calphad}
                        })
                else:
                    collection_calculation.update_one(
                        {'_id':ObjectId(task["_id"])}, 
                        {
                            '$set': {f"{source}.mixing_enthalpy": formation_enthalpy,
                                     f"{source}.mixing_enthalpy_for_calphad": formation_enthalpy_for_calphad}
                        })




