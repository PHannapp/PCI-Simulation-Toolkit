import itertools
import copy
import os
from src.config import n_hydrogen, n_POSCAR
from src.create_unblocked import create_unblocked_members, create_unblocked_members_solid
from MongoDB.connect import collection_calculation, collection_system, db

def create_reference_elements():
    systems = collection_system.find()
    for system in systems:
        if not system["elements"]:
            print("No elements in system!")
            return 
        for site in system["elements"]:
            if site == ["H,Va"]:
                site = ["H"]
            for element in site:
                if collection_calculation.find_one({"type": "reference", "elements": element}): # TODO maybe doesnt has to be from same system?
                    continue
                # Build the POSCAR path dynamically
                poscar_path = f"/home/pet69516/Dokumente/DFTforCALPHAD_home/POSCARS/{element.upper()}"
                
                # Check if POSCAR file exists at the specified path
                if os.path.exists(poscar_path):
                    with open(poscar_path, 'r') as poscar_file:
                        poscar_content = poscar_file.read()
                else:
                    string = f"No POSCAR for element {element}"
                    raise ValueError(string)

                task_data = {
                    # TODO automated determination of optimal DFT parameters (with ML?)
                    "system": system["system"],
                    "system_id": system["_id"],
                    "type": "reference",
                    "elements": element,
                    "POSCAR": poscar_content,
                    "cores": system["calculation_information"]["cores"],
                    "TO_RUN": False

                }
                # Inserting the document
                collection_calculation.insert_one(task_data)
                print(f"Inserted element {element} in reference_elements.")

def create_endmembers_new(system):        
    all_endmembers = []  # This will hold all the endmember combinations
    # Handle the special case where elements are combined like "H,Va"
    elements_lists = []
    for elements in system["elements"]:
        if elements == ["H,Va"]:
            elements = ['H']  # Split the combined elements
        elements_lists.append(elements)
    # Use itertools.product to generate all combinations (Cartesian product)
    combinations = list(itertools.product(*elements_lists))
    # Convert each tuple into a list of lists, i.e., [['Ce'], ['Al'], ['H'], ['H']] format
    combinations_as_lists = [[list([element]) for element in combo] for combo in combinations]
    # Append all combinations for this system to the main list
    all_endmembers.extend(combinations_as_lists)
    blocking_radius = 2.2
    for endmember in all_endmembers:
        if collection_calculation.find_one({"system_name": system["system"], "type": "endmember", "elements": endmember}): 
            continue
        for n_H in range(system["max_hydrogen"] + 1):
            poscars, blocking_radius = create_unblocked_members(system, endmember=endmember, n_H=n_H, n_POSCAR=n_POSCAR, blocking_radius=blocking_radius)
    
            for ith_poscar, poscar in enumerate(poscars):
                task_data = {
                    "system_name": system["system"],
                    "system_id": system["_id"],
                    "type": "endmember",
                    "elements": endmember,
                    "POSCAR": poscar,
                    "blocking_radius": blocking_radius,
                    "n_H": n_H,
                    "n_POSCAR": ith_poscar,
                    "cores": system["calculation_information"]["cores"],
                    "TO_RUN": False
                }
                # Inserting the document
                collection_calculation.insert_one(task_data)
            print(f"Inserted endmember {system['system']}H{n_H} in collection_endmembers.")


def create_mixing_new(system):
    all_mixing_species = []  # This will hold all the combinations with mixing species
    elements_lists = []
    for elements in system["elements"]:
        if elements == ["H,Va"]:
            elements = ['H','Va']  # Split the combined elements
        elements_lists.append(elements)
    # Identify sublattices with multiple elements
    for i, elements in enumerate(elements_lists):
        if elements == ['H', 'Va']: # already did as 'endmember'
            continue
        if len(elements) > 1:
            # Get all combinations of two elements
            combinations_of_two = list(itertools.combinations(elements, 2))
            for combination_of_two in combinations_of_two:
                elements_lists_copy = copy.deepcopy(elements_lists)
                del elements_lists_copy[i]
                combinations = list(itertools.product(*elements_lists_copy))
                for combination in combinations:
                    combination = [[element] for element in combination]
                    combination.insert(i, list(combination_of_two))
                    all_mixing_species.append(combination)
    for mixing_species in all_mixing_species:
        if collection_calculation.find_one({"system_name": system["system"], "type": "mixing", "elements": mixing_species}): 
            continue     
        # sort mixing sublattice alphabetically
        # create solid mixing with empty interstitial and with full interstitial
        poscarsii, metal_ratios, blocking_radii, n_H = create_unblocked_members_solid(system, endmember=mixing_species, n_POSCAR=n_POSCAR)
        if not poscarsii:
            print("Va endmember not yet calculated")
            continue
        for ith_ratio, metal_ratio in enumerate(metal_ratios):
            for ith_poscar, poscar in enumerate(poscarsii[ith_ratio]):
                task_data = {
                    "system_name": system["system"],
                    "system_id": system["_id"],
                    "type": "mixing",
                    "elements": mixing_species,
                    "POSCAR": poscar,
                    "blocking_radius": blocking_radii[ith_ratio],
                    "metal_ratio": metal_ratio,
                    "n_H": n_H,
                    "n_POSCAR": ith_poscar,
                    "cores": system["calculation_information"]["cores"],
                    "TO_RUN": False
                }
                # Inserting the document
                collection_calculation.insert_one(task_data)
                print(f"Inserted SQS {mixing_species} in collection_SQS.")