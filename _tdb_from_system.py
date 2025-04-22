from bson.objectid import ObjectId  # Required to work with _id
import os
from src.calculate import get_nested_value
from src.config import dS_H2
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution  # Import optimization functions
import numpy as np
import itertools
import csv
from MongoDB.connect import collection_calculation, collection_system, db

calculation_collection_name = collection_calculation.name

system_id = '678a29768c3e8eaa28454296'

def tdb_from_system(system_id=None, sources=None):
    if sources is None:
        sources = ["VASP"]
    for source in sources:
        if system_id is None:
            task_systems = collection_system.find({})
        else:
            task_systems = collection_system.find(
                {
                    '_id': ObjectId(system_id)
                }
            )
        for task_system in task_systems:
            database_name = ""
            tmp_str = ""
            unique_elements = []
            ## create the header string and a unique list of all different elements
            for s, sublattice in enumerate(task_system['elements']):
                for i, element in enumerate(sublattice):
                    if i != 0:
                        tmp_str += '-'
                    if element == 'H,Va':
                        tmp_str += 'H-Va'
                        unique_elements.append('H2')
                    else:
                        tmp_str += element
                        unique_elements.append(element)
                tmp_str += "_"
                tmp_str += str(task_system['multiplicity'][s])
                database_name += tmp_str
                tmp_str = "_"
            database_name += "_" + source.split('_')[0]
            unique_elements = list(set([element.upper() for element in unique_elements]))
            # 1. create file and name, import numpy
            final_string = "import numpy as np\n# reference -------------------------------------------------------------------------------\n"
            # 2. create ref
            for element in unique_elements:
                with open(os.path.join("CALPHAD_FILES", 'CONSTANTS', 'GHSER.py'), 'r') as file:
                    for line in file:
                        if f'GHSER{element}' in line:
                            # Append the line
                            final_string += line
                            # Try to get the next line, if it exists
                            next_line = next(file, None)
                            if next_line:
                                final_string += next_line
    
            final_string += "\n# Endmembers -------------------------------------------------------------------------------\n"
            # 3. create endmembers
            tasks = collection_calculation.find(
                {
                    'system_id': ObjectId(system_id),
                    'type': 'endmember',
                    f'{source}.formation_enthalpy_for_calphad': {'$exists': True}
                }
            )
    
            # Aggregation pipeline
            pipeline = [
                {
                    # Match documents based on the specified conditions
                    "$match": {
                        'system_id': ObjectId(system_id),
                        "type": "endmember",
                        f"{source}.formation_enthalpy_for_calphad": {"$exists": True}
                    }
                },
                {
                    "$group": {
                        "_id": {
                            "elements": "$elements",
                            "mixing_proportion": "$mixing_proportion"
                        },
                        "min_eff_energy": {"$min": f"${source}.formation_enthalpy_for_calphad"},
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
                        "from": calculation_collection_name,  # Performs a join between the current results and calculation_collection_name
                        "localField": "doc_id", # Matches the doc_id (stored in the previous step) with _id in calculation_collection_name.
                        "foreignField": "_id", # Stores the joined results in an array field "lowest_energy_doc"
                        "as": "lowest_energy_doc"
                    } # Purpose: Retrieves the full document from calculation_collection_name that corresponds to the minimum energy document.
                },
                {
                    "$unwind": { # Expands the "lowest_energy_doc" array so that we have a single document instead of an array.
                        "path": "$lowest_energy_doc",
                        "preserveNullAndEmptyArrays": True  # True ensures that if there’s no match in the lookup stage, the document is retained with null values.
                    } # Purpose: Ensures each document is flattened after the lookup.
                },
                {
                    "$replaceRoot": { # Replaces the current document structure with the full document from "lowest_energy_doc".
                        "newRoot": "$lowest_energy_doc"
                    } # Purpose: Returns the full document corresponding to the minimum enthalpy record instead of the intermediate aggregation results.
                }
            ]
    
            # Execute the pipeline
            results = list(collection_calculation.aggregate(pipeline))
            #def GCE_NI_V_V(T):
            #    return -197460 + (24.3084) * T + GHSERCE(T) + 5*GHSERNI(T)
            for task in results:
                line1 = 'def G'
                line2 = '   return '
                # for the calculation of the entropy
                amount_of_Hs = 0
                for sublattice, element in enumerate(task['elements']):
                    if sublattice != 0:
                        line1 += '_'
                    if element[0] == 'H':
                        line1 += 'H'
                        amount_of_Hs += task_system['multiplicity'][sublattice] / 2
                    elif element[0] == 'Va':
                        line1 += 'V'
                    else:
                        line1 += element[0].upper()
                        line2 += f'{task_system["multiplicity"][sublattice]}*GHSER{element[0].upper()}(T) + '
                line2 += f'{amount_of_Hs}*GHSERH2(T) + '
                hydrogen_entropy = amount_of_Hs * dS_H2
    
                source_parts = source.split(".")
                nested_source_value = get_nested_value(task, source_parts)
                formation_enthalpy = round(nested_source_value['formation_enthalpy_for_calphad'] * 1e3)
                line2 += f'({formation_enthalpy}) + ({hydrogen_entropy} * T)\n'
                line1 += '(T):\n'
                final_string += line1 + line2
    
            final_string += "\n# Interaction parameters -------------------------------------------------------------------\n"
    
            # 5. write interaction parameters in tdb
    
            tasks = collection_calculation.find(
                {
                    'system_id': ObjectId(system_id),
                    'type': 'mcsqs',
                    f'{source}.interaction_parameters': {'$exists': True}
                }
            )
            #def LCE_NI_HV_H_0(T):
            #    return -20632
            task_elements = []
            for task in tasks:
                source_parts = source.split(".")
                nested_source_value = get_nested_value(task, source_parts)
                if not nested_source_value.get("number_of_atoms_per_formula_unit"):
                    continue
                if task_elements != task['elements']:
                    task_elements = task['elements']
                else:
                    continue
                line1 = 'def L'
                # for the calculation of the entropy
                for sublattice, element in enumerate(task['elements']):
                    if sublattice != 0:
                        line1 += '_'
                    if element == ['H', 'Va']:
                        line1 += 'HV'
                    elif len(element) == 2:
                        line1 += element[0].upper()
                        line1 += element[1].upper()
                    elif element[0] == 'Va':
                        line1 += 'V'
                    else:
                        line1 += element[0].upper()
                for order, interaction_parameter in enumerate(nested_source_value["interaction_parameters"]):
                    line1_tmp = line1 + f'_{order}(T):\n'
                    line2 = f'   return {interaction_parameter }\n'#/ nested_source_value["number_of_atoms_per_formula_unit"]}\n'
                    final_string += line1_tmp + line2
            #def LCE_NI_HV_H_0(T):
            #    return -20632
    
    
            #### NOT YET
            # 6. calculate Entropy
    
            collection_system.update_one(
                {
                    '_id': ObjectId(system_id)
                },
                {
                    '$set': {f'{source}_tdb': final_string} 
                }
            )
            source_parts = source.split('.')
            source_str = '_'.join(source_parts)
            database_appendix = str(task_system['_id']) + '_' + source_str
            with open(os.path.join("CALPHAD_CALC", f'database_{database_appendix}.py'), 'w') as file:
                file.write(final_string)
            collection_system.update_one({'_id': task_system['_id']}, {'$set': {f"{source_str}_tdb": final_string}})
                
            print(final_string)

def create_pipeline(type, task_system, source, max_hydrogen, mixing=False):
    # Base match dictionary
    match_dict = {
        "$match": {
            "system_id": task_system["_id"],
            "type": type,
            f"{source}.formation_enthalpy_for_calphad": {"$exists": True},
            "n_H": {"$lt": max_hydrogen+1}
        }
    }
    # Conditionally add $or if type is "endmember"
    if type == "endmember" and not mixing:
        match_dict["$match"]["$or"] = [
            {"n_H": 0},
            {"n_H": max_hydrogen}
        ]
    pipeline = [
                match_dict,
                {
                    "$sort": {f"{source}.formation_enthalpy": 1}  # Sort by ascending energy
                },
                {
                    "$group": {
                        "_id": {
                            "elements": "$elements",
                            "n_H": "$n_H"
                        },
                        "min_eff_energy": {"$first": f"${source}.formation_enthalpy"},  # Select lowest enthalpy
                        "doc_id": {"$first": "$_id"}  # Capture ID of document with lowest energy
                    }
                },
                {
                    "$project": { # Extracts fields from _id (elements, mixing_proportion).
                        "elements": "$_id.elements",
                        "mixing_proportion": "$_id.mixing_proportion",
                        "min_eff_energy": 1, # Keeps min_eff_energy and doc_id.
                        "doc_id": 1
                    } # Purpose: Restructures the output for clarity before performing the next lookup.
                },
                {
                    "$lookup": {
                        "from": calculation_collection_name,  # Performs a join between the current results and calculation_collection_name
                        "localField": "doc_id", # Matches the doc_id (stored in the previous step) with _id in calculation_collection_name.
                        "foreignField": "_id", # Stores the joined results in an array field "lowest_energy_doc"
                        "as": "lowest_energy_doc"
                    } # Purpose: Retrieves the full document from calculation_collection_name that corresponds to the minimum energy document.
                },
                {
                    "$unwind": { # Expands the "lowest_energy_doc" array so that we have a single document instead of an array.
                        "path": "$lowest_energy_doc",
                        "preserveNullAndEmptyArrays": True  # True ensures that if there’s no match in the lookup stage, the document is retained with null values.
                    } # Purpose: Ensures each document is flattened after the lookup.
                },
                {
                    "$replaceRoot": { # Replaces the current document structure with the full document from "lowest_energy_doc".
                        "newRoot": "$lowest_energy_doc"
                    } # Purpose: Returns the full document corresponding to the minimum enthalpy record instead of the intermediate aggregation results.
                },
                {
                    "$group": {  # Grouping documents by the "elements" field
                        "_id": "$elements",
                        "grouped_docs": {"$push": "$$ROOT"}  # Collect all documents with the same elements into an array
                    }
                },
                {
                    "$project": {  # Optional: Rename the output structure
                        "_id": 0,  # Remove _id field
                        "elements": "$_id",
                        "documents": "$grouped_docs"
                    }
                }
            ]
    return pipeline

def sort_lists(first_list, second_list):
    sorted_pairs = sorted(zip(first_list, second_list))  # Sort based on first_list
    first_list_sorted, second_list_sorted = zip(*sorted_pairs)  # Unzip into separate lists
    return list(first_list_sorted), list(second_list_sorted)

def objective_function(fit, x_values, y_values):
    RMSE = 0
    for x, y in zip(x_values, y_values):
        fit_y = x*(1-x)*(fit[0]+fit[1]*(x-(1-x))+fit[2]*(x-(1-x))**2)
        if fit_y > y:
            RMSE += 10*(y - fit_y)**2
        else:
            RMSE += (y - fit_y)**2
    return (RMSE**(1/2)) / len(y_values)

def optimize_interaction(x_values, energies_trimmed):
    bounds = []
    bound = 1e7
    parameters = 3
    total = 3
    for i in range(parameters):
        bounds.append((-1*bound, bound))
    for i in range(total-parameters):
        bounds.append((0,0))
    result = differential_evolution(
        objective_function, args=(x_values, energies_trimmed), bounds=bounds, workers=1, disp=False, popsize=8, maxiter=4000
    )
    return result

def update_interaction_parameters(collection_system, task_system, new_elements, new_values):
    # Define the query to find the document with matching elements
    query = {
        '_id': task_system['_id'],
        'interaction_parameters.elements': new_elements  # Check if the elements already exist
    }

    # Define the update operation to overwrite the values if elements exist
    update_existing = { 
        '$set': {
            'interaction_parameters.$.values': list(new_values)  # Update the existing values
        }
    }

    # Try to update an existing entry
    result = collection_system.update_one(query, update_existing)

    # If no document was updated, append a new entry
    if result.matched_count == 0:
        new_entry = {
            "elements": new_elements,
            "values": list(new_values)
        }
        collection_system.update_one(
            {'_id': task_system['_id']},
            {'$push': {"interaction_parameters": new_entry}}
        )



def calculate_interaction_parameters_own(task_system, source, max_hydrogen):
    if source is None:
        source = "VASP"
    # search for H-Va
    interactions = ["endmember", "mixing"]
    for interaction in interactions:
        if interaction == "endmember":
            pipeline = create_pipeline("endmember", task_system, source, max_hydrogen, mixing=False)
        else:
            pipeline = create_pipeline("mixing", task_system, source, max_hydrogen)
        results = list(collection_calculation.aggregate(pipeline))
        for result in results:
            #if result["elements"][-1] == ["H"]:
            #    elements = result["elements"][:-1].append(["H,Va"])
            #else:
            #    elements = result["elements"]
            #mixing_proportion = []
            #energies = []
            #for task in result["documents"]:
            #    if result["elements"][-1] == ["H"]:
            #        mixing_proportion.append(task["n_H"] / task_system["max_hydrogen"])
            #    else:
            #        mixing_proportion.append(task["metal_ratio"])
            #    energies.append(task["VASP"]["formation_enthalpy_for_calphad"] * 1e3 * # per formula unti
            #                    task["number_of_atoms"] / task["VASP"]["number_of_atoms_per_formula_unit"]) # per supercell
            #mixing_proportion, energies = sort_lists(mixing_proportion, energies) # sort lists
            #energies_trimmed = [energies[i] - mixing_proportion[i]  * energies[-1] - (1 - mixing_proportion[i]) * energies[0] for i in range(len(energies))]
            #result = optimize_interaction(mixing_proportion, energies_trimmed)
            energies_dict = {"formation_enthalpy": {}}
            energies_dict_full = {"formation_enthalpy": {}}
            elements = result["elements"]
            tasks = collection_calculation.find(
                {
                    'system_name': task_system["system"],
                    "elements": elements,
                    "VASP.total_forces": {"$lt": 0.2},
                    "VASP.formation_enthalpy": {'$exists': True},
                    "n_H": {"$lt": max_hydrogen+1}
                }
            )
            for task in tasks:
                
                if interaction == "endmember":
                    n_H_atoms = n_H_atoms = task['n_H']
                    print(task['system_name'], task["n_H"], task["VASP"]["formation_enthalpy"])
                else:
                    n_H_atoms = 1 - task['metal_ratio'] # TODO
                    print(task["_id"], task['system_name'], 1 - task["metal_ratio"], task["VASP"]["formation_enthalpy_for_calphad"] * 1e3 * task["number_of_atoms"] / task["VASP"]["number_of_atoms_per_formula_unit"])
                if energies_dict_full["formation_enthalpy"].get(n_H_atoms):
                    energies_dict_full["formation_enthalpy"][n_H_atoms].append(task["VASP"]["formation_enthalpy_for_calphad"] * 1e3 * task["number_of_atoms"] / task["VASP"]["number_of_atoms_per_formula_unit"])
                else:
                    energies_dict_full["formation_enthalpy"][n_H_atoms] = [task["VASP"]["formation_enthalpy_for_calphad"] * 1e3 * task["number_of_atoms"] / task["VASP"]["number_of_atoms_per_formula_unit"]]
                if energies_dict["formation_enthalpy"].get(n_H_atoms):
                    if energies_dict["formation_enthalpy"][n_H_atoms] > task["VASP"]["formation_enthalpy_for_calphad"] * 1e3 * task["number_of_atoms"] / task["VASP"]["number_of_atoms_per_formula_unit"]:
                        energies_dict["formation_enthalpy"][n_H_atoms] = task["VASP"]["formation_enthalpy_for_calphad"] * 1e3 * task["number_of_atoms"] / task["VASP"]["number_of_atoms_per_formula_unit"]
                else:
                    energies_dict["formation_enthalpy"][n_H_atoms] = task["VASP"]["formation_enthalpy_for_calphad"] * 1e3 * task["number_of_atoms"] / task["VASP"]["number_of_atoms_per_formula_unit"]

            if interaction != "endmember":
                elements_Va = elements[:-1]
                elements_Va.append(["H"])
                combinations = list(itertools.product(*elements_Va))
                # Convert each tuple into a list of lists, i.e., [['Ce'], ['Al'], ['H'], ['H']] format
                endmembers = [[list([element]) for element in combo] for combo in combinations]
                count = 0
                for endmember in endmembers:
                    BREAK = False
                    task = collection_calculation.find_one(
                        {
                            'system_name': task_system["system"],
                            "elements": endmember,
                            "VASP.formation_enthalpy": {'$exists': True},
                            "n_H": 0
                        }
                    )
                    if not task:
                        print(f"Endmember {endmember}0 not found.")
                        BREAK = True
                        break
                    energies_dict["formation_enthalpy"][count] = task["VASP"]["formation_enthalpy_for_calphad"] * 1e3 * task["number_of_atoms"] / task["VASP"]["number_of_atoms_per_formula_unit"]
                    energies_dict_full["formation_enthalpy"][count] = [task["VASP"]["formation_enthalpy_for_calphad"] * 1e3 * task["number_of_atoms"] / task["VASP"]["number_of_atoms_per_formula_unit"]]
                    count += 1
                if BREAK:
                    break
            # Convert to sorted lists
            n_Hs, xs, energies = zip(*sorted(
                ((n_H, n_H, energy) for n_H, energy in energies_dict["formation_enthalpy"].items()),
                key=lambda x: x[0]  # Sort by x_H values
            ))
            if interaction == "endmember":
                energies_trimmed = [energies[i] - n_Hs[i] / n_Hs[-1]  * energies_dict["formation_enthalpy"][n_Hs[-1]] - (1 - n_Hs[i] / n_Hs[-1]) * energies_dict["formation_enthalpy"][0] for i in range(len(energies))]
            else:
                energies_trimmed = [energies[i] - n_Hs[i]  * energies_dict["formation_enthalpy"][n_Hs[-1]] - (1 - n_Hs[i]) * energies_dict["formation_enthalpy"][0] for i in range(len(energies))]
            energies_full_trimmed = []
            for n_H in n_Hs:
                if interaction == "endmember":
                    energy_trimmed_array = [energies_dict_full["formation_enthalpy"][n_H][i] - n_H / n_Hs[-1]  * energies_dict["formation_enthalpy"][n_Hs[-1]] - (1 - n_H / n_Hs[-1]) * energies_dict["formation_enthalpy"][0] for i in range(len(energies_dict_full["formation_enthalpy"][n_H]))]
                else:
                    energy_trimmed_array = [energies_dict_full["formation_enthalpy"][n_H][i] - n_H  * energies_dict["formation_enthalpy"][n_Hs[-1]] - (1 - n_H) * energies_dict["formation_enthalpy"][0] for i in range(len(energies_dict_full["formation_enthalpy"][n_H]))]
                energies_full_trimmed.append(energy_trimmed_array)

            if interaction == "endmember":
                x_values = [n_H / max(n_Hs) for n_H in n_Hs]
            else:
                x_values = [n_H for n_H in n_Hs]
            for ith_nH, n_H in enumerate(n_Hs):
                if interaction == "endmember":
                    n_H = n_H / max(n_Hs)
                plt.plot([n_H for _ in energies_full_trimmed[ith_nH]], energies_full_trimmed[ith_nH], 'o', color='grey')    
            plt.plot(x_values, energies_trimmed, color='green', label="Minimum energy trajectory")
            # Fitting equation: =x*(1-x)*(fit[0]+fit[1]*(x-(1-x))+fit[2]*(x-(1-x))**2)
            opt_result = optimize_interaction(x_values, energies_trimmed)
            y_fit_values = [x*(1-x)*(opt_result.x[0]+opt_result.x[1]*(x-(1-x))+opt_result.x[2]*(x-(1-x))**2) for x in np.linspace(0,1,100)]
            label_appendix = ""
            for res in opt_result.x[:3]:
                label_appendix += f"{int(res)},  "
            plt.plot(np.linspace(0,1,100), y_fit_values, color='red', label=f"Polynomial fit: {label_appendix}")                       
            plt.xlabel("Hydrogen occupancy")
            plt.ylabel("J / mol (POSCAR)")
            plt.legend()
            if result["elements"][-1] == ["H"]:
                elements = result["elements"][:-1]
                elements.append(["H,Va"])
            else:
                elements = result["elements"]
            el_str = "-"
            for sl in elements:
                for el in sl:
                    el_str += el
                el_str += "-"
            with open(os.path.join(r"/home/pet69516/Dokumente/DFTforCALPHAD_home/PLOTS/", "csvs", f"{task_system['system']}_{el_str}_max{str(max_hydrogen)}.csv"), mode='w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                # Write the first row (database module names)
                writer.writerow(["All"])
                for ith_nH, n_H in enumerate(n_Hs):
                    if interaction == "endmember":
                        n_H = n_H / max(n_Hs)
                    for x, p in zip([n_H for _ in energies_full_trimmed[ith_nH]], energies_full_trimmed[ith_nH]):
                        writer.writerow([x, p])

                writer.writerow(["Min"])
                for x, p in zip(x_values, energies_trimmed):
                    writer.writerow([x, p])

                writer.writerow(opt_result.x)
                for x, p in zip(np.linspace(0,1,100), y_fit_values):
                    writer.writerow([x, p])
            plt.savefig(r"/home/pet69516/Dokumente/DFTforCALPHAD_home/PLOTS/" + task_system["system"] + el_str + "_max" + str(max_hydrogen) + ".png")
            plt.close()
            update_interaction_parameters(collection_system, task_system, elements, opt_result.x)
        
def tdb_from_system_own(system_id=None, sources=None, max_hydrogen=None):
    if sources is None:
        sources = ["VASP"]
    for source in sources:
        if system_id is None:
            task_systems = collection_system.find({'max_hydrogen': {'$exists': True}})#, "interaction_parameters": {'$exists': False}})
        else:
            task_systems = collection_system.find(
                {
                    '_id': ObjectId(system_id),
                    'max_hydrogen': {'$exists': True}
                }
            )
        for task_system in task_systems:
            if max_hydrogen is None:
                max_hydrogen = task_system["max_hydrogen"]
            unique_elements = []
            ## create the header string and a unique list of all different elements
            for s, sublattice in enumerate(task_system['elements']):
                for i, element in enumerate(sublattice):
                    if element == 'H,Va':
                        unique_elements.append('H2')
                    else:
                        unique_elements.append(element)
            unique_elements = list(set([element.upper() for element in unique_elements]))
            # 1. create file and name, import numpy
            final_string = "import numpy as np\n# reference -------------------------------------------------------------------------------\n"
            # 2. create ref
            for element in unique_elements:
                with open(os.path.join("CALPHAD_FILES", 'CONSTANTS', 'GHSER.py'), 'r') as file:
                    for line in file:
                        if f'GHSER{element}' in line:
                            # Append the line
                            final_string += line
                            # Try to get the next line, if it exists
                            next_line = next(file, None)
                            if next_line:
                                final_string += next_line
    
            final_string += "\n# Endmembers -------------------------------------------------------------------------------\n"
            # 3. create endmembers
            # Aggregation pipeline
            pipeline = create_pipeline("endmember", task_system, source, max_hydrogen)
            # Execute the pipeline
            results = list(collection_calculation.aggregate(pipeline))
            #def GCE_NI_V_V(T):
            #    return -197460 + (24.3084) * T + GHSERCE(T) + 5*GHSERNI(T)
            multiplicity = None
            for result in results:
                for task in result["documents"]:
                    print(task["elements"])
                    print(task["n_H"])
                    # determine the multiplicity of the supercell
                    if multiplicity is None:
                        multiplicity = task_system["structure"]["simple"]["multiplicities"] # is the unit cell multiplicity
                        n_metal_atoms = task["number_of_atoms"] - task["n_H"]
                        multiplier = n_metal_atoms / sum(multiplicity) # determine number of unitcells
                        multiplicity = [m * multiplier for m in multiplicity]
                        multiplicity.append(max_hydrogen)
                    line1 = 'def G'
                    line2 = '   return '
                    # for the calculation of the entropy
                    amount_of_Hs = 0
                    for sublattice, element in enumerate(task['elements']):
                        if sublattice != 0:
                            line1 += '_'
                        if element[0] == 'H' and task["n_H"] == 0:
                            line1 += 'V'
                        elif element[0] == 'H' and task["n_H"] == max_hydrogen:
                            line1 += 'H'
                            amount_of_Hs += multiplicity[sublattice] / 2
                        elif element[0] == 'Va':
                            line1 += 'V'
                        else:
                            line1 += element[0].upper()
                            line2 += f'{multiplicity[sublattice]}*GHSER{element[0].upper()}(T) + '
                    line2 += f'{amount_of_Hs}*GHSERH2(T) + '
                    hydrogen_entropy = amount_of_Hs * dS_H2

                    source_parts = source.split(".")
                    nested_source_value = get_nested_value(task, source_parts)
                    formation_enthalpy = round(nested_source_value['formation_enthalpy_for_calphad'] * 1e3) # per task["VASP"]["number_of_atoms_per_formula_unit"]
                    formation_enthalpy *= task["number_of_atoms"] / task["VASP"]["number_of_atoms_per_formula_unit"]

                    line2 += f'({formation_enthalpy}) + ({hydrogen_entropy} * T)\n'
                    line1 += '(T):\n'
                    final_string += line1 + line2
    
            final_string += "\n# Interaction parameters -------------------------------------------------------------------\n"
    
            # 5. write interaction parameters in tdb
            calculate_interaction_parameters_own(task_system, source, max_hydrogen)
            task_system = collection_system.find_one({"_id": task_system["_id"]})
            #def LCE_NI_HV_H_0(T):
            #    return -20632
            for interactions in task_system["interaction_parameters"]:
                interaction_elements = interactions["elements"]
                interaction_parameters = interactions["values"]
                line1 = 'def L'
                # for the calculation of the entropy
                for sublattice, element in enumerate(interaction_elements):
                    if sublattice != 0:
                        line1 += '_'
                    if element == ['H', 'Va']:
                        line1 += 'HV'
                    elif element == ['H,Va']:
                        line1 += 'HV'
                    elif len(element) == 2:
                        line1 += element[0].upper()
                        line1 += element[1].upper()
                    elif element[0] == 'Va':
                        line1 += 'V'
                    else:
                        line1 += element[0].upper()
                for order, interaction_parameter in enumerate(interaction_parameters):
                    line1_tmp = line1 + f'_{order}(T):\n'
                    line2 = f'   return {interaction_parameter }\n'#/ nested_source_value["number_of_atoms_per_formula_unit"]}\n'
                    final_string += line1_tmp + line2
            #def LCE_NI_HV_H_0(T):
            #    return -20632
    
    
            #### NOT YET
            # 6. calculate Entropy
    
            source_parts = source.split('.')
            source_str = '_'.join(source_parts)
            database_appendix = str(task_system["system"]) + '_' + source_str + '_max' + str(max_hydrogen)
            with open(os.path.join("CALPHAD_CALC", f'database_{database_appendix}.py'), 'w') as file:
                file.write(final_string)
            collection_system.update_one({'_id': task_system['_id']}, {'$set': {f"{source_str}_max{str(max_hydrogen)}_tdb": final_string}})
                
            print(final_string)

for n_H in [20,32]:
    print(n_H)
    tdb_from_system_own('67ae084088ff415532bc6d57', max_hydrogen=n_H)