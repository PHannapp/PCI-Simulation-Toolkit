from pymongo import MongoClient
from bson.objectid import ObjectId  # Required to work with _id
import matplotlib.pyplot as plt
from pymatgen.core import Structure
import py3Dmol
from calculate import fit_interaction_parameters, fit_interaction_parameters2
import re

# Connect to MongoDB
client = MongoClient('mongodb://localhost:27017/')
db = client.DFTforCALPHAD
collection_system = db.System
collection_calculation = db.calculation

calculation_collection_name = collection_calculation.name

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

import os 
task_dir = r"/home/pet69516/Dokumente/DFTforCALPHAD/VASPResults"
def print_energies(system, keyword1, keyword2, name = None, energy='energy'):
    tasks = collection_calculation.find(
        {
            'system': system,
            keyword1: {'$exists': True}
        }
    )
    for task in tasks:
        if name:
            if name not in task['name']:
                continue

        source_parts = keyword1.split(".")
        nested_source_value = get_nested_value(task, source_parts)
        if nested_source_value.get(keyword2):
            ##print(f"{task['name']}: {task[keyword1][keyword2]['energy']}")
            #print(f"{task['name']}: {task[keyword1][keyword2]['formation_enthalpy']}")
            if nested_source_value[keyword2].get(energy):
                if task.get('name'):
                    print(task['name'],f"{nested_source_value[keyword2][energy]}")
                else:
                    print(task['_id'],",",task["elements"],",",f"{nested_source_value[keyword2][energy]}")
                
                #if task[keyword1][keyword2].get("CONTCAR"):
                #    with open(os.path.join(task_dir, f"{task['name']}.poscar"), 'w') as POTCARfile:
                #        POTCARfile.write(task[keyword1][keyword2]['CONTCAR'])
                #else:
                #    print("NO_CONTCAR!")
            else:
                if task.get('name'):
                    print(task['name'],"ERROR")
                else:
                    print(task['_id'],"ERROR")

def print_energies_from_OSZICAR(system, keyword1, keyword2):
    tasks = collection_calculation.find(
        {
            'system': system,
            keyword1: {'$exists': True}
        }
    )
    for task in tasks:
        if r"C14" not in task['name']:
            continue
        if task[keyword1].get(keyword2):
            ##print(f"{task['name']}: {task[keyword1][keyword2]['energy']}")
            #print(f"{task['name']}: {task[keyword1][keyword2]['formation_enthalpy']}")
            energy = None
            if task[keyword1][keyword2].get('OSZICAR_EFF'):
                oszicar_array = task[keyword1][keyword2]['OSZICAR_EFF'].splitlines()
                for line in oszicar_array:
                    if "E0=" in line:
                        line_parts = line.split()
                        for m, line_part in enumerate(line_parts):
                            if line_part == "E0=":
                                energy = float(line_parts[m + 1]) / task['number_of_atoms']
                                break
                print(task['name'], energy)

def plot_endmember_energies(keyword1, keyword2, force, exclude, save=None):
    system_tasks = collection_system.find({})
    for system_task in system_tasks:
        tasks = collection_calculation.find(
            {
                'system_id': system_task['_id'],
                'type': 'endmember',
                '$or': [
                    {f'{keyword1}.total_forces': {'$lt': force}}, 
                    {f'{keyword1}.total_forces': None}
                ],
                f'{keyword1}.{keyword2}': {'$exists': True}
            }
        )

        # Aggregation pipeline
        pipeline = [
            {
                # Match documents based on the specified conditions
                "$match": {
                    "system_id": system_task['_id'],
                    "type": "endmember",
                    '$or': [
                        {f'{keyword1}.total_forces': {'$lt': force}}, 
                        {f'{keyword1}.total_forces': None}
                    ],
                    f"{keyword1}.{keyword2}": {"$exists": True}
                }
            },
            {
                "$group": {
                    "_id": {
                        "elements": "$elements",
                        "mixing_proportion": "$mixing_proportion"
                    },
                    "min_eff_energy": {"$min": f"${keyword1}.{keyword2}"},
                    "doc_id": {"$first": "$_id"}  # Capture the ID of one document with the lowest energy
                }
            },
            {
                # Debug: Add a stage to output intermediate results
                "$project": {
                    "elements": "$_id.elements",
                    "mixing_proportion": "$_id.mixing_proportion",
                    "min_eff_energy": 1,
                    "doc_id": 1
                }
            },
            {
                "$lookup": {
                    "from": calculation_collection_name,  # Replace with the actual collection name
                    "localField": "doc_id",
                    "foreignField": "_id",
                    "as": "lowest_energy_doc"
                }
            },
            {
                "$unwind": {
                    "path": "$lowest_energy_doc",
                    "preserveNullAndEmptyArrays": True  # Set to True to identify mismatches
                }
            },
            {
                "$replaceRoot": {
                    "newRoot": "$lowest_energy_doc"
                }
            }
        ]

        # Execute the pipeline
        results = list(collection_calculation.aggregate(pipeline))

        multiplicity = system_task['multiplicity']
        result_dict = {}
        for task in results:
            BREAK = False
            key = ""
            x = 0
            for i, element in enumerate(task['elements']):
                if element[0] == exclude:
                    BREAK = True
                    break
                if element[0] == "H":
                    x += multiplicity[i]
                elif element[0] != 'Va':
                    if multiplicity[i] == 1:
                        key += element[0]
                    else:
                        key += element[0] + str(multiplicity[i])
            if BREAK:
                continue
            # Prepare the data to add to the result dictionary
            source_parts = keyword1.split(".")
            nested_source_value = get_nested_value(task, source_parts)
            new_entry = [x, nested_source_value[keyword2]]  # [x, energy or another value]
            # If the key already exists, append the new entry to the list
            if key in result_dict:
                result_dict[key].append(new_entry)
            else:
                # If the key doesn't exist, create a new list with the new entry
                result_dict[key] = [new_entry]
        # Now plot the results using matplotlib
        if len(result_dict) == 0:
            return
        plt.figure(figsize=(10, 6))  # Define figure size
        for key, values in result_dict.items():
            # Sort the values based on x (the first element in the sublist)
            values.sort(key=lambda pair: pair[0])

            # Extract x and y values from the sorted arrays
            x_values = [item[0] for item in values]
            y_values = [item[1] for item in values]

            # Plot each key's data with a label
            plt.plot(x_values, y_values, marker='o', linestyle='-', label=f"{key}-H")

        # Set plot labels and title
        plt.xlabel('H-atoms', fontsize=14)
        plt.ylabel(f'{keyword2} ({keyword1})', fontsize=14)
        plt.title(f'Endmember Energies for System {system_task["system"]}', fontsize=16)

        # Add legend
        plt.legend(loc='best')
        if save is not None:
            # Replace problematic characters in the system name
            safe_system_name = re.sub(r'[<>:"/\\|?*]', '_', system_task['system'])
            output_filename = f"{safe_system_name}_{keyword1}_{keyword2}_endmember.png"  # Define your file name
            plt.savefig(f"{save}/{output_filename}", dpi=300)  # Save plot with 300 dpi
        else:
            # Display the plot
            plt.show()


def plot_mixing_energies(system, keyword1, keyword2, force, exclude=None, save=None):
    tasks = collection_calculation.find(
        {
            'system': system,
            'type': 'mcsqs',
            '$or': [
                {f'{keyword1}.total_forces': {'$lt': force}}, 
                {f'{keyword1}.total_forces': None}
            ],
            f'{keyword1}.{keyword2}': {'$exists': True}
        }
    )

    # Aggregation pipeline
    pipeline = [
        {
            # Match documents based on the specified conditions
            "$match": {
                "system": system,
                "type": "mcsqs",
                '$or': [
                    {f'{keyword1}.total_forces': {'$lt': force}}, 
                    {f'{keyword1}.total_forces': None}
                ],
                f"{keyword1}.{keyword2}": {"$exists": True}
            }
        },
        {
            "$group": {
                "_id": {
                    "elements": "$elements",
                    "mixing_proportion": "$mixing_proportion"
                },
                "min_eff_energy": {"$min": f"${keyword1}.{keyword2}"},
                "doc_id": {"$first": "$_id"}  # Capture the ID of one document with the lowest energy
            }
        },
        {
            # Debug: Add a stage to output intermediate results
            "$project": {
                "elements": "$_id.elements",
                "mixing_proportion": "$_id.mixing_proportion",
                "min_eff_energy": 1,
                "doc_id": 1
            }
        },
        {
            "$lookup": {
                "from": calculation_collection_name,  # Replace with the actual collection name
                "localField": "doc_id",
                "foreignField": "_id",
                "as": "lowest_energy_doc"
            }
        },
        {
            "$unwind": {
                "path": "$lowest_energy_doc",
                "preserveNullAndEmptyArrays": True  # Set to True to identify mismatches
            }
        },
        {
            "$replaceRoot": {
                "newRoot": "$lowest_energy_doc"
            }
        }
    ]

    # Execute the pipeline
    results = list(collection_calculation.aggregate(pipeline))

    if not results:
        print("No tasks retrieved.")
        return
    system_task = collection_system.find_one(
        {
            'system': system
        })
    if not system_task:
        print("No systems retrieved.")
        return
    multiplicity = system_task['multiplicity']
    # Initialize an array of dictionaries for each sublattice
    result_dict_array = [{} for _ in multiplicity]
    for task in results:
        BREAK = False
        key = ""
        x = 0
        for i, sublattice in enumerate(task['elements']):
            if exclude in sublattice:
                BREAK = True
                break
            if len(sublattice) == 2:
                mixing_sublattice = i
                if multiplicity[i] == 1:
                    key += f"({sublattice[0]}ₓ{sublattice[1]}₁₋ₓ)"
                else:
                    key += f"({sublattice[0]}ₓ{sublattice[1]}₁₋ₓ)" + str(multiplicity[i])
            else:
                if multiplicity[i] == 1:
                    key += sublattice[0]
                else:
                    key += sublattice[0] + str(multiplicity[i])
        if BREAK:
            continue
        

        # Prepare the data to add to the result dictionary
        source_parts = keyword1.split(".")
        nested_source_value = get_nested_value(task, source_parts)
        new_entry = [task['mixing_proportion'], nested_source_value[keyword2]]  # [x, energy or another value]
        # If the key already exists, append the new entry to the list
        if key in result_dict_array[mixing_sublattice]:
            result_dict_array[mixing_sublattice][key].append(new_entry)
        else:
            # If the key doesn't exist, create a new list with the new entry
            result_dict_array[mixing_sublattice][key] = [new_entry]
    for sublattice, dict in enumerate(result_dict_array):
        if len(dict) == 0:
            continue
        for key in dict:
            # Append [0, 0] and [1, 0] to the existing list for each key in the sublattice dictionary
            dict[key].append([0, 0])
            dict[key].append([1, 0])
        # Now plot the results using matplotlib
        plt.figure(figsize=(10, 6))  # Define figure size
        for key, values in dict.items():
            # Sort the values based on x (the first element in the sublist)
            values.sort(key=lambda pair: pair[0])

            # Extract x and y values from the sorted arrays
            x_values = [item[0] for item in values]
            y_values = [item[1] for item in values]

            # Plot each key's data with a label
            plt.plot(x_values, y_values, marker='o', linestyle='-', label=f"{key}")

        # Set plot labels and title
        plt.xlabel(f'x in AₓB₁₋ₓ', fontsize=14)
        plt.ylabel(f'{keyword2} ({keyword1})', fontsize=14)
        plt.title(f'Mixing Energies for mixing in sublattice {sublattice} in System {system}', fontsize=16)

        # Add legend
        plt.legend(loc='best')
        if save is not None:
            output_filename = f"{keyword1}_{keyword2}_mixing_sublattice_{sublattice}.png"  # Define your file name
            plt.savefig(f"{save}/{output_filename}", dpi=300)  # Save plot with 300 dpi
        else:
            # Display the plot
            plt.show()

def fit_mixing_energies_to_interaction_parameters(system, keyword1, bound, force, parameter_space, save=None):
    tasks = collection_calculation.find(
        {
            'system': system,
            'type': 'mcsqs',
            '$or': [
                {f'{keyword1}.total_forces': {'$lt': force}}, 
                {f'{keyword1}.total_forces': None}
            ],
            f'{keyword1}.mixing_enthalpy_for_calphad': {'$exists': True}
        }
    )
    system_task = collection_system.find_one(
        {
            'system': system
        })
    multiplicity = system_task['multiplicity']
    # Initialize an array of dictionaries for each sublattice
    result_dict_array = [{} for _ in multiplicity]
    key_dict_array = [{} for _ in multiplicity]
    # check for multiple nth

    # Aggregation pipeline
    pipeline = [
        {
            # Match documents based on the specified conditions
            "$match": {
                "system": system,
                "type": "mcsqs",
                '$or': [
                    {f'{keyword1}.total_forces': {'$lt': force}}, 
                    {f'{keyword1}.total_forces': None}
                ],
                f"{keyword1}.mixing_enthalpy_for_calphad": {"$exists": True}
            }
        },
        {
            "$group": {
                "_id": {
                    "elements": "$elements",
                    "mixing_proportion": "$mixing_proportion"
                },
                "min_eff_energy": {"$min": f"${keyword1}.mixing_enthalpy_for_calphad"},
                "doc_id": {"$first": "$_id"}  # Capture the ID of one document with the lowest energy
            }
        },
        {
            # Debug: Add a stage to output intermediate results
            "$project": {
                "elements": "$_id.elements",
                "mixing_proportion": "$_id.mixing_proportion",
                "min_eff_energy": 1,
                "doc_id": 1
            }
        },
        {
            "$lookup": {
                "from": calculation_collection_name,  # Replace with the actual collection name
                "localField": "doc_id",
                "foreignField": "_id",
                "as": "lowest_energy_doc"
            }
        },
        {
            "$unwind": {
                "path": "$lowest_energy_doc",
                "preserveNullAndEmptyArrays": True  # Set to True to identify mismatches
            }
        },
        {
            "$replaceRoot": {
                "newRoot": "$lowest_energy_doc"
            }
        }
    ]

    # Execute the pipeline
    results = list(collection_calculation.aggregate(pipeline))






    for task in results:
        key = str(task['elements'])
        for i, sublattice in enumerate(task['elements']):
            if len(sublattice) == 2:
                mixing_sublattice = i
                ISHV = True
                if sublattice != ['H', 'Va']:
                    ISHV = False
        # Prepare the data to add to the result dictionary
        source_parts = keyword1.split(".")
        nested_source_value = get_nested_value(task, source_parts)
        new_entry = [task['mixing_proportion'], nested_source_value['mixing_enthalpy_for_calphad']]  # [x, energy or another value]
        # If the key already exists, append the new entry to the list
        if key in result_dict_array[mixing_sublattice]:
            result_dict_array[mixing_sublattice][key].append(new_entry)        
        else:
            # If the key doesn't exist, create a new list with the new entry
            result_dict_array[mixing_sublattice][key] = [new_entry]
        if key not in key_dict_array[mixing_sublattice]:
            key_dict_array[mixing_sublattice][key] = task['elements']
    for sublattice, dict in enumerate(result_dict_array):
        if len(dict) == 0:
            continue
        #for key in dict:
        #    # Append [0, 0] and [1, 0] to the existing list for each key in the sublattice dictionary
        #    dict[key].append([0, 0])
        #    dict[key].append([1, 0])
        if save is not None:
            plt.figure(figsize=(10, 6))  # Define figure size
        for key, values in dict.items():
            # Sort the values based on x (the first element in the sublist)
            values.sort(key=lambda pair: pair[0])            # Extract x and y values from the sorted arrays
            x_values = [item[0] for item in values]
            y_values = [item[1]*1e3 for item in values]
            if False:#ISHV:
                interaction_parameters, x_fit_values, y_fit_values= fit_interaction_parameters(x_values, y_values, key, bound, parameter_space, save)
                y_values = [y - interaction_parameters[3] for y in y_values]
            else:
                interaction_parameters, x_fit_values, y_fit_values= fit_interaction_parameters2(x_values, y_values, key, bound, parameter_space, save)
            interaction_parameters = [int(i) for i in interaction_parameters]
            # Apply the update to the document
            collection_calculation.update_many(
                {
                    'system_id': system_task['_id'],
                    'elements': key_dict_array[sublattice][key]
                },
                {
                    '$set': {f'{keyword1}.interaction_parameters': interaction_parameters[:3]} 
                }
            )    
            if save is not None:
                # Plot each key's data with a label
                # Plot the data and retrieve the color automatically
                line, = plt.plot(x_values, y_values, marker='o')

                # Use the same color for the fit line
                plt.plot(x_fit_values, y_fit_values, color=line.get_color(), label=f"{key}")

        if save is not None:
            # Set plot labels and title
            plt.xlabel(f'x in AₓB₁₋ₓ', fontsize=14)
            plt.ylabel(f'{keyword1}-mixing_enthalpy_for_calphad)', fontsize=14)
            plt.title(f'Mixing Energies for mixing in sublattice {sublattice} in System {system}', fontsize=16)
            plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
            # Add legend
            plt.tight_layout()
            output_filename = f"{keyword1}_interaction_parameter_{sublattice}.png"  # Define your file name
            plt.savefig(f"{save}/{output_filename}", dpi=300)  # Save plot with 300 dpi

def fit_mixing_energies_to_interaction_parameters2(system, keyword1, bound, force, parameter_space, save=None):
    tasks = collection_calculation.find(
        {
            'system': system,
            'type': 'mcsqs',
            '$or': [
                {f'{keyword1}.total_forces': {'$lt': force}}, 
                {f'{keyword1}.total_forces': None}
            ],
            f'{keyword1}.mixing_enthalpy_for_calphad': {'$exists': True}
        }
    )
    system_task = collection_system.find_one(
        {
            'system': system
        })
    multiplicity = system_task['multiplicity']
    # Initialize an array of dictionaries for each sublattice
    result_dict_array = [{} for _ in multiplicity]
    key_dict_array = [{} for _ in multiplicity]
    for task in tasks:
        key = str(task['elements'])
        for i, sublattice in enumerate(task['elements']):
            if len(sublattice) == 2:
                mixing_sublattice = i
        # Prepare the data to add to the result dictionary
        source_parts = keyword1.split(".")
        nested_source_value = get_nested_value(task, source_parts)
        new_entry = [task['mixing_proportion'], nested_source_value['mixing_enthalpy_for_calphad']]  # [x, energy or another value]
        # If the key already exists, append the new entry to the list
        if key in result_dict_array[mixing_sublattice]:
            result_dict_array[mixing_sublattice][key].append(new_entry)        
        else:
            # If the key doesn't exist, create a new list with the new entry
            result_dict_array[mixing_sublattice][key] = [new_entry]
        if key not in key_dict_array[mixing_sublattice]:
            key_dict_array[mixing_sublattice][key] = task['elements']
    for sublattice, dict in enumerate(result_dict_array):
        if len(dict) == 0:
            continue
        #for key in dict:
        #    # Append [0, 0] and [1, 0] to the existing list for each key in the sublattice dictionary
        #    dict[key].append([0, 0])
        #    dict[key].append([1, 0])
        if save is not None:
            plt.figure(figsize=(10, 6))  # Define figure size
        for key, values in dict.items():
            # Sort the values based on x (the first element in the sublist)
            values.sort(key=lambda pair: pair[0])            # Extract x and y values from the sorted arrays
            x_values = [item[0] for item in values]
            y_values = [item[1]*1e3 for item in values]

            interaction_parameters, x_fit_values, y_fit_values= fit_interaction_parameters(x_values, y_values, key, bound, parameter_space, save)
            interaction_parameters = [int(i) for i in interaction_parameters]
            # Apply the update to the document
            collection_calculation.update_many(
                {
                    'system_id': system_task['_id'],
                    'elements': key_dict_array[sublattice][key]
                },
                {
                    '$set': {f'{keyword1}.interaction_parameters': interaction_parameters} 
                }
            )
            if save is not None:
                # Plot each key's data with a label
                # Plot the data and retrieve the color automatically
                line, = plt.plot(x_values, y_values, marker='o', label=f"{key}")

                # Use the same color for the fit line
                plt.plot(x_fit_values, y_fit_values, color=line.get_color(), label=f"{key}")

        if save is not None:
            # Set plot labels and title
            plt.xlabel(f'x in AₓB₁₋ₓ', fontsize=14)
            plt.ylabel(f'{keyword1}-mixing_enthalpy_for_calphad)', fontsize=14)
            plt.title(f'Mixing Energies for mixing in sublattice {sublattice} in System {system}', fontsize=16)

            # Add legend
            plt.legend(loc='best')
            output_filename = f"{keyword1}_interaction_parameter_{sublattice}.png"  # Define your file name
            plt.savefig(f"{save}/{output_filename}", dpi=300)  # Save plot with 300 dpi

def show_strcuture():
    tasks = collection_calculation.find(
        {
            "mixing_proportion": 0.25,
            "elements.1": ["Al", "Ni"],
            "elements.2": ["Va"],
            "elements.3": ["H"]
        })
    for task in tasks:
        # Load the POSCAR file using pymatgen
        structure = Structure.from_str(task["CHGNet_results"]["final_structure"], fmt='poscar')
        # Export the structure to CIF format (since XYZ is not supported directly)
        cif_str = structure.to(fmt="cif")
        # Print structure to check if it's loaded properly
        print(structure)
        # Print CIF string to debug
        print(cif_str)
        # Create a py3Dmol viewer
        view = py3Dmol.view(width=800, height=600)
        view.addModel(cif_str, "cif")  # Use CIF format for py3Dmol
        view.setStyle({"stick": {}})
        view.zoomTo()
        
        # Show the structure in the browser
        view.show()