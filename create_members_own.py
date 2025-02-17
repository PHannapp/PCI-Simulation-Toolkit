import copy
import time
import numpy as np
import itertools
import hashlib


class Atom:
    def __init__(self, element, position, label):
        self.element = element
        self.position = position
        self.label = label
        self.n_nearest_neighbors = None  # Initialize the attribute

class Structure:
    def __init__(self, lattice_vectors, atoms):
        self.lattice_vectors = lattice_vectors
        self.atoms = atoms

        vac_atoms = []
        metal_atoms = []
        for atom in atoms:
            if "Vac" in atom.label:
                vac_atoms.append(atom)
            else:
                metal_atoms.append(atom)

        self.vac_atoms = vac_atoms
        self.metal_atoms = metal_atoms

# Helper function to compute a hash of an array or sub-array
def compute_hash(array):
    array_string = str(array)  # Convert array or sub-array to string
    return hashlib.md5(array_string.encode()).hexdigest()  # Use md5 hash

# Function to sort arrays by their hash values
def sort_by_hash(arrays):
    return sorted(arrays, key=lambda arr: compute_hash(arr))

# Function to compare two arrays of arrays
def are_arrays_equal(array1, array2):
    # Sort both arrays by their hash values
    sorted_array1 = sort_by_hash(array1)
    sorted_array2 = sort_by_hash(array2)
    
    # Compare the sorted arrays
    return np.array_equal(sorted_array1, sorted_array2)

# Function to compare two dictionaries, ignoring the integer in the key tuple
def compare_dicts(dict1, dict2, tolerance=1e-6):
    # Group by the atom labels (ignoring the integer)
    grouped_dict1 = group_by_label(dict1)
    grouped_dict2 = group_by_label(dict2)
    #print("\n\n GOUPED DICT 1###############################")
    #print("grouped_dict1", grouped_dict1)
    #print("\n\n GOUPED DICT 2###############################")
    #print("grouped_dict2", grouped_dict2)
    # Compare grouped dictionaries
    if set(grouped_dict1.keys()) != set(grouped_dict2.keys()):
        return False  # If the labels (keys without the integers) are not the same
    
    # Compare the arrays associated with each label
    for label in grouped_dict1:
        arrays1 = sort_by_hash(grouped_dict1[label])
        arrays2 = sort_by_hash(grouped_dict2[label])
        #print("\n\n###############################")
        #print(arrays1)
        #print("\n\n###############################")
        #print(arrays2)
        # Compare each pair of arrays element-wise
        for (distances1, atoms1), (distances2, atoms2) in zip(arrays1, arrays2):
            # Compare the distance arrays with a tolerance
            if not np.allclose(distances1, distances2, atol=tolerance):
                return False  # If the distances don't match
            
            # Compare the atom label arrays
            if not np.array_equal(atoms1, atoms2):
                return False  # If the atom labels don't match
    
    return True

# Group the dictionary by the label part of the key (ignoring the integer)
def group_by_label(label_distances_dict):
    grouped_dict = {}
    #print("#############")
    #input(label_distances_dict)
    for (label, _), (distances, atoms) in label_distances_dict.items():
        if label not in grouped_dict:
            grouped_dict[label] = []
        grouped_dict[label].append([distances, atoms])
    return grouped_dict

def compute_distances_with_periodic(sites, lattice_vectors, supercell_size, site_type):
    """
    Compute the minimum distances between atoms in a supercell, considering periodic boundary conditions.
    """
    vacant_sites_positions = sites[site_type]["positions"]
    number_of_atoms = len(vacant_sites_positions)
    distances = np.zeros((number_of_atoms, number_of_atoms))
    sites[site_type]["shift_vectors"] = np.empty((number_of_atoms, number_of_atoms, 3))
    # Define all possible translations (shifts in any direction by -1, 0, or 1 unit cell)
    shifts = list(itertools.product([-1, 0, 1], repeat=3))
    
    for i in range(number_of_atoms):
        for j in range(i + 1, number_of_atoms):  # Only calculate for i < j to exploit symmetry
            # Calculate the minimum distance by considering the wrap-around effect in all directions
            min_distance = np.inf
            for shift in shifts:                
                # Convert shift and supercell_size to arrays and multiply element-wise
                scaled_shift = np.array(shift) * np.array(supercell_size)
                # Apply shift to the second atom's position
                shifted_site_j = vacant_sites_positions[j] + np.dot(scaled_shift, lattice_vectors)
                # Calculate the Euclidean distance between the atom i and the shifted atom j
                distance = np.linalg.norm(vacant_sites_positions[i] - shifted_site_j)
                # Keep track of the minimum distance found
                if distance < min_distance:
                    min_distance = distance
                    # TODO MAYBE IS SHIFT UNNECESSARY TO STORE
                    sites[site_type]["shift_vectors"][i][j] = shift
                    sites[site_type]["shift_vectors"][j][i] = -np.array(shift)  # Symmetric shift for j,i
            distances[i, j] = min_distance
            distances[j, i] = min_distance  # Copy the value to the symmetric position
    
    # Set distances that are exactly 0 (i.e., diagonal elements) to infinity
    np.fill_diagonal(distances, np.inf)
    sites[site_type]["distances"] = distances  
    return sites

def setup_sites(supercell_size, atoms, lattice_vectors):
    number_of_unit_cells = np.prod(supercell_size)
    total_number_of_atoms = number_of_unit_cells * len(atoms)
    vacant_sites = np.zeros((total_number_of_atoms, 3)) # = [[x, y, z], [x, y, z], ...]
    vacant_sites_labels = np.empty((total_number_of_atoms), dtype=object) # = [[x, y, z], [x, y, z], ...]
    vacant_translation_vectors = np.zeros((total_number_of_atoms, 3)) # = [[x, y, z], [x, y, z], ...]
    # initialize all sites, by placing the atoms in the original and translated unit cells and directly multiply by the lattice vector# Define individual ranges for each axis
    x_max, y_max, z_max = supercell_size
    x_range = np.arange(0, x_max) 
    y_range = np.arange(0, y_max) 
    z_range = np.arange(0, z_max) 
    # Generate the translation vectors using itertools.product
    translation_vectors = np.array(list(itertools.product(x_range, y_range, z_range)))
    count = 0
    for translation_vector in translation_vectors: 
        for ith_atom_in_unitcell, vac_atom in enumerate(atoms):
            vacant_sites[count] = np.dot((vac_atom.position + translation_vector), lattice_vectors)
            vacant_translation_vectors[count] = translation_vector
            vacant_sites_labels[count] = vac_atom.label
            count += 1
    #distances = np.linalg.norm(vacant_sites[:, np.newaxis] - vacant_sites, axis=2)
    sites = {"vacant": 
               {"positions": vacant_sites, 
                "label": vacant_sites_labels,
                "translation_vectors": vacant_translation_vectors, 
                "shift_vectors": None, 
                "distances": None, 
                "idx": np.arange(0, total_number_of_atoms) },
             "occupied": 
               {"positions": np.empty((0, 3)), 
                "label": np.empty((0)), 
                "translation_vectors": np.empty((0, 3)), 
                "shift_vectors": np.empty((0, 3)), 
                "distances": None, 
                "idx": np.array(()) },
             "blocked":
               {"positions": np.empty((0, 3)), 
                "label": np.empty((0)), 
                "translation_vectors": np.empty((0, 3)), 
                "shift_vectors": np.empty((0, 3)),
                "distances": None, 
                "idx": np.array(()) } }
    sites = compute_distances_with_periodic(sites, lattice_vectors, supercell_size, "vacant")
    return sites, total_number_of_atoms

def eliminate_duplicates(occupied_positions, occupied_labels, sites_array, tolerance=1e-6):
    """
    Eliminate duplicate configurations from occupied_positions by comparing each position directly.
    """
    unique_positions = [occupied_positions[0]]  # To store unique configurations
    unique_labels = [occupied_labels[0]]  # To store unique configurations
    unique_sites_array = [sites_array[0]]  # To store unique configurations

    # Loop over remaining configurations starting from the second one
    for i, position_array_to_test in enumerate(occupied_positions[1:]):
        # we first think not all are the same -> its unique:
        ARE_SAME = False
        # Loop over the already stored unique positions: if one whole position_array_to_test_against has same positions as position_array_to_test: break and NOT append
        for position_array_to_test_against in unique_positions:
            result = are_arrays_equal(position_array_to_test_against, position_array_to_test)
            if result:
                #(f"Alarm, arrays are the same: {position_array_to_test_against}, {position_array_to_test}")
                ARE_SAME = True
                break
        if not ARE_SAME:
            unique_positions.append(position_array_to_test)
            unique_labels.append(occupied_labels[i + 1])
            unique_sites_array.append(sites_array[i + 1])
            #print(1)
            #input(position_array_to_test)
        #else:
            #print(2)
            #input(position_array_to_test)

    return np.array(unique_positions), np.array(unique_labels), np.array(unique_sites_array)

def  eliminate_symmetry_duplicates(occupied_positions, occupied_labels, unique_sites_array, total_number_of_metal_atoms, tolerance=1e-6):
    """
    Eliminate duplicate configurations from occupied_positions by comparing each position directly.
    """
    # 1. i have distances of every atom to every other
    # 2. i have everies atom wyckoff position
    # 3. i take only the hydrogen atoms distances
    # 4. i sort the distances, to have the minimum 4 (or 6 maybe)
    # 5. two conditions must be met: 
    #       1. all distance arrays have to be the same
    #       2. if two distance arrays are the same, the elements that point to the distances must be the same
    unique_positions = [occupied_positions[0]]  # To store unique configurations
    # TODO create a label: distances-matrix dict - array (matrix is len(sorted_occupied_positions[i] x len(sorted_occupied_positions[i] size))
    label_distances_dict_array = np.empty((len(occupied_positions)), dtype=object)
    for indx in range(len(occupied_positions)):
        label_distances_dict = {}
        # 3. i take only the hydrogen atoms distances
        positions = unique_sites_array[indx]["occupied"]["positions"][total_number_of_metal_atoms:]
        distances = unique_sites_array[indx]["occupied"]["distances"][total_number_of_metal_atoms:]
        labels = unique_sites_array[indx]["occupied"]["label"][total_number_of_metal_atoms:]
        for i, (label, position) in enumerate(zip(labels, positions)):
            label_distances_dict[(label, i)] = [distances[i], unique_sites_array[indx]["occupied"]["label"]]
        label_distances_dict_array[indx] = label_distances_dict
    #print(label_distances_dict_array)
    for label_distances_dict in label_distances_dict_array:
        for key, value in label_distances_dict.items():
            distances, other_atoms = value
            # Sort the distances from minimum to maximum and reorder other_atoms accordingly
            sorted_indices = np.argsort(distances)  # Get indices that would sort the distances
            sorted_distances = np.array(distances)[sorted_indices]  # Apply sorting to distances
            sorted_other_atoms = np.array(other_atoms)[sorted_indices]  # Apply sorting to atom names
            
            # Replace the dictionary values with the sorted ones
            label_distances_dict[key] = [sorted_distances, sorted_other_atoms]

    unique_label_distances_dict_array = [label_distances_dict_array[0]] 
    # Loop over remaining configurations starting from the second one
    for i, label_distances_dict_to_test in enumerate(label_distances_dict_array[1:], start=1):
        #print("------------------------------- NEW TO TEST")
        IS_UNIQUE = True  # Assume this configuration is unique

        # Loop over all previously stored unique configurations
        for label_distances_dict_to_test_against in unique_label_distances_dict_array:
            if compare_dicts(label_distances_dict_to_test, label_distances_dict_to_test_against, tolerance):
                IS_UNIQUE = False  # If we find a match, mark as not unique
                break
        #print("+++++++++++++++++++++++++ ALL AGAINST TESTED")
        if IS_UNIQUE:
            unique_positions.append(occupied_positions[i])
            unique_label_distances_dict_array.append(label_distances_dict_to_test)

    return np.array(unique_positions)

def setup_analysis(supercell_size, structure):
    occupied_ratios = []
    unique_labels = list(set(atom.label for atom in structure.vac_atoms))
    total_number_of_atoms_per_label = {}
    occupation_ratios_per_labels = {}
    for label in unique_labels:
        occupation_ratios_per_labels[label] = []
        total_number_of_atoms_per_label[label] = 0
    for vac_atom in structure.vac_atoms:
        total_number_of_atoms_per_label[vac_atom.label] += 1
    for label in unique_labels:
        total_number_of_atoms_per_label[label] *= np.prod(supercell_size)
    return occupied_ratios, total_number_of_atoms_per_label, occupation_ratios_per_labels

def occupation_per_label(sites, structure):
    # Extract the occupation in each label:
    unique_labels = list(set(atom.label for atom in structure.vac_atoms))
    occupation_labels = {}
    for label in unique_labels:
        occupation_labels[label] = []
    
    n_vac_atoms = len(structure.vac_atoms)
    for idx, position in zip(sites["occupied"]["idx"], sites["occupied"]["positions"]):
        occupation_labels[structure.vac_atoms[int(idx) % n_vac_atoms].label].append(position)
    return occupation_labels

def delete_and_append_for_matrix(sites, i_blocked, i_occupied):
    # Append the ith position and idx to the specified category
    for i in i_blocked:
        if i == i_occupied:
            where_to_append = "occupied"
        else:
            where_to_append = "blocked"
        sites[where_to_append]["positions"] = np.vstack(
            [sites[where_to_append]["positions"], sites["vacant"]["positions"][i]]
        )
        sites[where_to_append]["shift_vectors"] = np.vstack(
            [sites[where_to_append]["shift_vectors"], sites["vacant"]["shift_vectors"][i]]
        )
        sites[where_to_append]["idx"] = np.append(sites[where_to_append]["idx"], sites["vacant"]["idx"][i])
        sites[where_to_append]["label"] = np.append(sites[where_to_append]["label"], sites["vacant"]["label"][i])
        #print("Position: ", sites["vacant"]["positions"][i], " \nLabel: ", sites["vacant"]["label"][i])
        #input()
    # Delete the ith element from 'vacant'
    sites["vacant"]["positions"] = np.delete(sites["vacant"]["positions"], i_blocked, axis=0)
    sites["vacant"]["idx"] = np.delete(sites["vacant"]["idx"], i_blocked)
    sites["vacant"]["label"] = np.delete(sites["vacant"]["label"], i_blocked)
    sites["vacant"]["shift_vectors"] = np.delete(sites["vacant"]["shift_vectors"], i_blocked, axis=0)
    
    sites["vacant"]["distances"] = np.delete(sites["vacant"]["distances"], i_blocked, axis=0)
    sites["vacant"]["distances"] = np.delete(sites["vacant"]["distances"], i_blocked, axis=1)

def fill_cell_with_ratio_hydrogen(sites, blocking_radius, vac_atoms, total_number_of_atoms_per_label, ratio_dict, number_of_atoms):
    n_vacant_sites = len(sites["vacant"]["positions"])
    n_sites_to_reach_ratio = number_of_atoms
    t_now = time.time()
    iteration = 0
    n_occupied = 0
    total_ratio_per_label = {}
    n_vac_atoms = len(vac_atoms)
    for label, value in total_number_of_atoms_per_label.items():
        total_ratio_per_label[label] = 0
    next_label = None
    penalty = 0
    while n_occupied < n_sites_to_reach_ratio:
        if n_vacant_sites <= 0:
            raise ValueError("Error, ratio could not be reached")
        #if iteration == 10:
            #print(f"Time for 10 / {n_vacant_sites} possible iterations: ", time.time() - t_now)
        iteration += 1
        # choose ith atom of vacant sites (already excluded blocked sites)
        ith_atom = np.random.randint(0, n_vacant_sites)
        # check if the ith site is from the correct label
        idx = sites["vacant"]["idx"][ith_atom]
        label = vac_atoms[int(idx) % n_vac_atoms].label
        if next_label is None or  next_label == label:
            penalty = 0
            total_ratio_per_label[label] += 1 / total_number_of_atoms_per_label[label]
            # calculate what the next label should be, by determining the farest ratio
            differences = {key: (ratio_dict[key] - total_ratio_per_label[key]) for key in total_ratio_per_label}
            next_label = max(differences, key = differences.get)
        else:
            if penalty >= 100000:
                return None, False
            penalty += 1
            continue
        # Compute the Euclidean distances from the ith_atom to all other sites
        indices = np.where(sites["vacant"]["distances"][ith_atom] < blocking_radius)[0]
        #distances = np.linalg.norm(sites["vacant"]["positions"] - sites["vacant"]["positions"][ith_atom], axis=1)
        # Get the indices where the distance is less than the threshold
        #indices = np.where(distances < blocking_radius)[0]
        delete_and_append_for_matrix(sites, np.concatenate(([ith_atom], indices)), ith_atom)
        n_occupied += 1
        n_vacant_sites -= (len(indices) + 1)
    return sites, True

def create_endmember_with_defined_occupation_ratio(supercell_size, blocking_radius, number_of_iterations, ratio_dict, structure):
    occupied_ratios, total_number_of_atoms_per_label, occupation_ratios_per_labels = setup_analysis(supercell_size, structure)
    _, total_number_of_atoms_per_label, _ = setup_analysis(supercell_size, structure)
    sites, total_number_of_atoms =  setup_sites(supercell_size, structure.vac_atoms, structure.lattice_vectors)
    number_of_atoms = sum([ratio_dict[key] * total_number_of_atoms_per_label[key] for key in ratio_dict.keys()])
    copied_sites_array = np.empty((number_of_iterations), dtype=object)
    if round(number_of_atoms, 4) % 1 != 0:
        #return []
        raise ValueError(f"Ratio * total_number_of_Vac_atoms {total_number_of_atoms} is not possible, leading to non-int {number_of_atoms}.")
    for label in total_number_of_atoms_per_label:
        if round(total_number_of_atoms_per_label[label] * ratio_dict[label], 4) % 1 != 0:
            #return []
            raise ValueError(f"For label {label} ratio {ratio_dict[label]} * total_number_of_Vac_atoms {total_number_of_atoms_per_label[label]} is not possible, leading to non-int {total_number_of_atoms_per_label[label] * ratio_dict[label]}.")
    number_of_atoms = int(round(number_of_atoms))
    occupied_positions = np.empty((number_of_iterations, number_of_atoms, 3))
    occupied_labels = np.empty((number_of_iterations, number_of_atoms), dtype=object)
    iteration = 0  # Initialize iteration counter
    while iteration < number_of_iterations:
        print("Iteration: ", iteration, " / ", int(number_of_iterations))
        try:
            copied_sites, SUCCESS = fill_cell_with_ratio_hydrogen(copy.deepcopy(sites), blocking_radius, structure.vac_atoms, total_number_of_atoms_per_label, ratio_dict, number_of_atoms)
            if not SUCCESS:
                print("Could not find solution!")
                occupied_positions = np.delete(occupied_positions, -1, axis=0)  # Removes the last row along the first axis
                occupied_labels = np.delete(occupied_labels, -1, axis=0)  # Removes the last row along the first axis
                copied_sites_array = np.delete(copied_sites_array, -1, axis=0)  # Removes the last row along the first axis
                number_of_iterations -= 1  # Decrease the total number of iterations
                continue
            occupation_labels = occupation_per_label(copied_sites, structure)
            for label, array in occupation_labels.items():
                occupation_ratios_per_labels[label].append(len(array) / total_number_of_atoms_per_label[label] * 100)
        except ValueError as e:
            print(e)
            occupied_positions = np.delete(occupied_positions, -1, axis=0)  # Removes the last row along the first axis
            occupied_labels = np.delete(occupied_labels, -1, axis=0)  # Removes the last row along the first axis
            copied_sites_array = np.delete(copied_sites_array, -1, axis=0)  # Removes the last row along the first axis
            number_of_iterations -= 1  # Decrease the total number of iterations
            continue
        # Do not increment iteration, so it retries the same iteration
        occupied_positions[iteration] = copied_sites["occupied"]["positions"]
        occupied_labels[iteration] = copied_sites["occupied"]["label"]
        copied_sites_array[iteration] = copied_sites
        iteration += 1  # Only increment if successful
    return occupied_positions, occupied_labels, copied_sites_array

def extract_initialize(system, endmember, ratio):
    lattice_vectors = np.array(system['structure']["lattice_parameters"])
    ratio_dict = {}
    atoms = []
    vac_count = 0
    for sublattice, positions in enumerate(system['structure']['atomic_positions']):
        if endmember[sublattice] == ["Va"]:
            vac_count += 1
            continue
        elif endmember[sublattice] == ["H"]:
            label = f"Vac{sublattice}"
            ratio_dict[label] = ratio[vac_count]
            vac_count += 1
        else:
            label = f"{endmember[sublattice]}{sublattice}"
        element = endmember[sublattice]
        for position in positions:
            #pos = np.dot(np.array(position), lattice_vectors)
            pos = position
            atom = Atom(element, pos, label)
            atoms.append(atom)
    structure = Structure(lattice_vectors, atoms)
    return structure, ratio_dict    

def write_poscars(system, structure, metal_sites, symmetry_unique_occupied_positions, endmember, ratio, mixing_ratio):
    # Save optimized radii to a CSV file
    n_unitcell = np.prod(np.array(system['supercell_size']))
    for axis in range(len(system['supercell_size'])):
        for ith_vector, vector in enumerate(structure.lattice_vectors[axis]):
            structure.lattice_vectors[axis][ith_vector] = vector * system['supercell_size'][axis]
    
    name = ""
    mixing_sublattice = None
    for member in endmember:
        if len(member) == 2:
            name += member[0]
            name += member[1]
            continue
        name += member[0]
    basearray = [name, "1.0"]
    for vector in structure.lattice_vectors:
        basearray.append(' '.join(map(str, vector)))
    name_correct = ""
    for sublattice, member in enumerate(endmember):
        if member == ["Va"]:
            continue
        if len(member) == 2:
            ### prepare for positions
            mixing_sublattice = sublattice
            length_metal_sublattices = []
            for ith_endmember, metal_endmember in enumerate(endmember):
                if metal_endmember in [["H"], ["Va"], [["H"], ["Va"]]]:
                    break
                length_metal_sublattices.append(len(system['structure']['atomic_positions'][ith_endmember]))
            number_of_first_mixing_metal_atoms = length_metal_sublattices[mixing_sublattice] * n_unitcell * mixing_ratio

            name_correct += member[0] + " "
            name_correct += member[1] + " "
            continue
        name_correct += member[0] + " "
        if member == ["H"]:
            break
    basearray.append(name_correct)
    numbers = ""
    vac_sl = 0
    vac_count = 0
    for sublattice, positions in enumerate(system['structure']['atomic_positions']):
        if endmember[sublattice] == ["Va"]:
            vac_sl += 1
            continue
        elif endmember[sublattice] == ["H"]:
            vac_count += int(round(len(positions) * ratio[vac_sl]))
            vac_sl += 1
        else:
            if len(endmember[sublattice]) == 2:
                numbers += str(int(round(len(positions) * mixing_ratio * n_unitcell))) + " "
                numbers += str(int(round(len(positions) * (1 - mixing_ratio) * n_unitcell))) + " "
                continue
            numbers += str(int(round(len(positions) * n_unitcell))) + " "
    if vac_count:
        numbers += str(int(round(vac_count * n_unitcell))) + " "
    basearray.append(numbers)
    basearray.append('Cartesian')

    basestring = '\n'.join(map(str, basearray)) + '\n'

    number_of_metal_atoms_total = len(metal_sites["vacant"]["positions"])
    if number_of_metal_atoms_total / n_unitcell % 1 != 0:
        #return []
        raise ValueError(f"Something is wrong with number {number_of_metal_atoms_total} / {n_unitcell} leading to non-int {number_of_metal_atoms_total / n_unitcell}")
    number_of_metal_atoms_per_formula_unit = int(round(number_of_metal_atoms_total / n_unitcell))
    poscar_strings = []
    for enum, positions in enumerate(symmetry_unique_occupied_positions):
        string = basestring
        if mixing_sublattice is not None:
            positions_per_element_metal = [[] for _ in range(len(length_metal_sublattices))]
            for unit_cell in range(n_unitcell):
                nth_element = 0
                for sublattice, n_elements_in_sublattice in enumerate(length_metal_sublattices):
                    for _ in range(n_elements_in_sublattice):
                        positions_per_element_metal[sublattice].append(metal_sites["vacant"]["positions"][unit_cell * number_of_metal_atoms_per_formula_unit + nth_element])
                        nth_element += 1
                    
            # random decision of mixing metal elements
            mixing_positions = [[], []]
            number_of_second_mixing_metal_atoms = (len(positions_per_element_metal[mixing_sublattice]) - number_of_first_mixing_metal_atoms) 
            for ith_atom in range(len(positions_per_element_metal[mixing_sublattice])):
                # if all the rest positions should go in the second mixing element
                if len(positions_per_element_metal[mixing_sublattice]) - ith_atom == number_of_first_mixing_metal_atoms - len(mixing_positions[0]):
                    mixing_positions[1].append(positions_per_element_metal[mixing_sublattice][ith_atom])
                # if all the rest positions should go in the first mixing element
                elif len(positions_per_element_metal[mixing_sublattice]) - ith_atom == number_of_second_mixing_metal_atoms - len(mixing_positions[1]):
                    mixing_positions[0].append(positions_per_element_metal[mixing_sublattice][ith_atom])
                    # randomly choose
                else:
                    random = np.random.random()
                    if random > mixing_ratio:
                        mixing_positions[1].append(positions_per_element_metal[mixing_sublattice][ith_atom])
                    else:
                        mixing_positions[0].append(positions_per_element_metal[mixing_sublattice][ith_atom])
            # replace the mixiung positons
            if mixing_sublattice == len(positions_per_element_metal) - 1:
                positions_per_element_metal = positions_per_element_metal[:mixing_sublattice] + mixing_positions
            else:
                positions_per_element_metal = positions_per_element_metal[:mixing_sublattice] + mixing_positions + positions_per_element_metal[mixing_sublattice+1:]
            
            for sublattice in positions_per_element_metal:
                for position in sublattice:
                    string += ' '.join(map(str, position)) + "\n"
        else:
            for unit_cell in range(n_unitcell):
                for position in range(number_of_metal_atoms_per_formula_unit):
                    string += ' '.join(map(str, metal_sites["vacant"]["positions"][unit_cell * number_of_metal_atoms_per_formula_unit + position])) + "\n"
        for position in positions:
            string += ' '.join(map(str, position)) + "\n"
        poscar_strings.append(string)
    return poscar_strings

def create_unblocked_members(system, endmember, mixing_ratio=None):
    # get the .NN file
    blocking_radius = 2.2
    iterations = 10
    ratio = system['ratios']
    if mixing_ratio:
        # mixing member
        vac_site = 0
        for site in range(len(endmember)):
            if endmember[site] == ["H", "Va"]:
                ratio[vac_site] *= mixing_ratio
                endmember[site] = ["H"]
                break
            elif endmember[site] in [["H"], ["Va"]]:
                vac_site += 1
    supercell_size = system['supercell_size']
    structure, ratio_dict = extract_initialize(system, endmember, ratio)
    saved_occupied_positions = []
    saved_occupied_labels = []
    saved_sites_array = []
    while True:
        occupied_positions, occupied_labels, sites_array = create_endmember_with_defined_occupation_ratio(supercell_size, blocking_radius, iterations, ratio_dict, structure)
        saved_occupied_positions.extend(occupied_positions)
        saved_occupied_labels.extend(occupied_labels)
        saved_sites_array.extend(sites_array)
        if not len(saved_occupied_positions):
            blocking_radius -= 0.05
            print(f"################################################################# \n\nOnly {len(occupied_positions)} results where found. New blocking radius: ", blocking_radius, "\n\n\n")
        elif len(saved_occupied_positions) < iterations:
            print(f"################################################################# \n\nOnly {len(occupied_positions)} results where found. Repeat.\n\n\n")
        else:
            break
    t0 = time.time()
    unique_occupied_positions, unique_occupied_labels, unique_sites_array = eliminate_duplicates(saved_occupied_positions[:iterations], saved_occupied_labels[:iterations], saved_sites_array[:iterations])
    t1 = time.time()
    print(f"{len(unique_occupied_positions)} out of  {iterations} results are unique in {t1 - t0} seconds.\n")
    metal_sites, total_number_of_metal_atoms =  setup_sites(supercell_size, structure.metal_atoms, structure.lattice_vectors)
    for unique_sites in unique_sites_array:
        for key in metal_sites["vacant"].keys():
            if key not in ["distances", "shift_vectors", "idx"]:
                unique_sites["occupied"][key] = np.concatenate((metal_sites["vacant"][key], unique_sites["occupied"][key]))
        unique_sites = compute_distances_with_periodic(unique_sites, structure.lattice_vectors, supercell_size, "occupied")
    t2 = time.time()
    print(f"compute_distances_with_periodic done in {t2 - t1} seconds. Start eliminating symmetry duplicates.")
    symmetry_unique_occupied_positions = eliminate_symmetry_duplicates(unique_occupied_positions, unique_occupied_labels, unique_sites_array, total_number_of_metal_atoms)
    t3 = time.time()
    print(f"eliminate_symmetry_duplicates done in {t3 - t2} seconds.\n")
    #write_csv(output_dir, unique_occupied_positions, metal_sites, NN_file_path, "6h16h2_1supercell_111", ratio, supercell_size)
    print(f"{len(symmetry_unique_occupied_positions)} out of  {len(occupied_positions)} results are symmetrical unique.\n")

    poscars = write_poscars(system, structure, metal_sites, symmetry_unique_occupied_positions, endmember, ratio, mixing_ratio)
    return poscars, blocking_radius