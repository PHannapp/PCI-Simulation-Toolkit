from src.create_database_entries_from_system import create_reference_elements, create_endmembers_new, create_mixing_new
import os
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from src.make_Wyckoff_db import read_wyckoff, get_structure
from src.VoidSize import void_size
from src.read import change_poscar
from pymatgen.io.vasp import Poscar
from MongoDB.connect import collection_calculation, collection_system, db

# Save the current working directory
original_directory = os.getcwd()
directory = os.path.join(original_directory, "POSCARS", 'SystemPoscars')
name = 'La2Fe22Ti2'

def insertSystem(directory, name):
    system = collection_system.find_one({'system': name})
    if system:
        raise ValueError("Name already in system.")

    poscar_path = os.path.join(directory, name)
    with open(poscar_path, 'r') as poscar_file:
        poscar_string = poscar_file.read()

    poscar = Poscar.from_file(poscar_path)  # Load POSCAR
    structure = poscar.structure  # Extract Structure object
    # Get the space group
    sga = SpacegroupAnalyzer(structure)
    symmetry = sga.get_symmetry_dataset()
    structure_simple = get_structure(structure, symmetry["std_types"]) # input([1,2,2,2,2,2]), output: {positions:[[[0,0,0]], [[1,1,1],[1,1,1],[1,1,1],[1,1,1],[1,1,1]]], multiplicity:[1,4]}
    structure_complex = get_structure(structure, symmetry["equivalent_atoms"]) # [1,2,2,3,3,3], output: {positions:[[[0,0,0]], [[1,1,1],[1,1,1]], [[2,2,2],[2,2,2],[2,2,2]]], multiplicity:[1,2,3]}
    wyckoff = read_wyckoff("src/Wyckoff.csv")
    wyckoff_positions = wyckoff[symmetry.hall_number-1]["wyckoff"]
    change_poscar(poscar_path, symmetry["wyckoffs"])
    ### 1. Read metal poscar
    #create_NN()
    NN_string = void_size(poscar_path, wyckoff_positions)
    
    task_data = {
        "system": name,
        "elements": [["LA"], ["FE"], ["TI"], ["H,Va"]],
        "poscar": poscar_string,
        "NN": NN_string,
        "structure": {"simple":
            {"atomic_positions": structure_simple["positions"],
            "multiplicities": structure_simple["multiplicities"]
            },
            "complex":
            {"atomic_positions": structure_complex["positions"],
            "multiplicities": structure_complex["multiplicities"]
            }},
        "lattice_parameters": [list(vector) for vector in structure.lattice.matrix],
        "calculation_information":
        {
            "cutoff": [400],
            "k_point_density": [0.22] , 
            "cores": 24  
            },
        "n_metal_atoms": 12,
        "max_hydrogen": 32
    }
    # Inserting the document
    collection_system.insert_one(task_data)

def add_element(name=None):
    if name:
        systems = collection_system.find({'system': name})
    else:
        systems = collection_system.find({})
    for system in systems:
        # check if its a new system
        if not system.get("max_hydrogen"):
            continue
        ## check if element is in accordance with already stored elements    
        #if system["elements"] and len(element) != len(system["elements"][0]):
        #    raise ValueError("Structure of new element is not in accordance with existing elements")
        #if element in system["elements"]:
        #    raise ValueError("Element already in system.")
        #system_elements = copy.deepcopy(system['elements'])
        #system_elements.append(element)
        #collection_system.update_one({'_id': system['_id']}, {'$set': {'elements': system_elements}})
        #collection_system.update_one({'_id': system['_id']}, {'$push': {'elements': element}})
        create_reference_elements()
        create_endmembers_new(system)
        # check for reference element
        #for sl in element:
        #    task = collection_calculation.find_one({'elements': sl[0], 'type': 'reference'})
        #    if not task:
        #        create_reference_elements()

        # TODO check for mixing!
        create_mixing_new(system)


#insertSystem(directory, name)
add_element() # element has to be array of sublattices, ALSO if 

