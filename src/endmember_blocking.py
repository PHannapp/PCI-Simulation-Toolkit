import numpy as np
import os
from src.read import read_poscar
from src.calculate_NN import seperate_vac_atoms
import itertools
import copy
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import time

subfolder_name = "STEP3_NN"
colors = list(mcolors.TABLEAU_COLORS.values())  # List of pre-defined colors from Matplotlib

def compute_distances_with_periodic(sites, lattice_vectors, supercell_size):
    """
    Compute the minimum distances between atoms in a supercell, considering periodic boundary conditions.
    """
    vacant_sites_positions = sites["vacant"]["positions"]
    number_of_atoms = len(vacant_sites_positions)
    distances = np.zeros((number_of_atoms, number_of_atoms))
    
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
                    sites["vacant"]["shift_vectors"][i][j] = shift
                    sites["vacant"]["shift_vectors"][j][i] = -np.array(shift)  # Symmetric shift for j,i
            distances[i, j] = min_distance
            distances[j, i] = min_distance  # Copy the value to the symmetric position
    
    # Set distances that are exactly 0 (i.e., diagonal elements) to infinity
    np.fill_diagonal(distances, np.inf)
    sites["vacant"]["distances"] = distances
    
    return sites

def setup_sites(supercell_size, vac_atoms, lattice_vectors):
    number_of_unit_cells = np.prod(supercell_size)
    total_number_of_atoms = number_of_unit_cells * len(vac_atoms)
    vacant_sites = np.zeros((total_number_of_atoms, 3)) # = [[x, y, z], [x, y, z], ...]
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
        for ith_atom_in_unitcell, vac_atom in enumerate(vac_atoms):
            vacant_sites[count] = np.dot((vac_atom.position + translation_vector), lattice_vectors)
            vacant_translation_vectors[count] = translation_vector
            count += 1
    #distances = np.linalg.norm(vacant_sites[:, np.newaxis] - vacant_sites, axis=2)
    sites = {"vacant": 
               {"positions": vacant_sites, 
                "translation_vectors": vacant_translation_vectors, 
                "shift_vectors": np.empty((total_number_of_atoms, total_number_of_atoms, 3)), 
                "distances": None, 
                "idx": np.arange(0, total_number_of_atoms) },
             "occupied": 
               {"positions": np.empty((0, 3)), 
                "translation_vectors": np.empty((0, 3)), 
                "idx": np.array(()) },
             "blocked":
               {"positions": np.empty((0, 3)), 
                "translation_vectors": np.empty((0, 3)), 
                "idx": np.array(()) } }
    sites = compute_distances_with_periodic(sites, lattice_vectors, supercell_size)
    return sites, total_number_of_atoms

def setup_analysis(vac_atoms, supercell_size):
    occupied_ratios = []
    unique_labels = list(set(atom.label for atom in vac_atoms))
    total_number_of_atoms_per_label = {}
    occupation_ratios_per_labels = {}
    for label in unique_labels:
        occupation_ratios_per_labels[label] = []
        total_number_of_atoms_per_label[label] = 0
    for vac_atom in vac_atoms:
        total_number_of_atoms_per_label[vac_atom.label] += 1
    for label in unique_labels:
        total_number_of_atoms_per_label[label] *= np.prod(supercell_size)
    return occupied_ratios, total_number_of_atoms_per_label, occupation_ratios_per_labels

def occupation_per_label(sites, vac_atoms):
    # Extract the occupation in each label:
    unique_labels = list(set(atom.label for atom in vac_atoms))
    occupation_labels = {}
    for label in unique_labels:
        occupation_labels[label] = []
    
    n_vac_atoms = len(vac_atoms)
    for idx, position in zip(sites["occupied"]["idx"], sites["occupied"]["positions"]):
        occupation_labels[vac_atoms[int(idx) % n_vac_atoms].label].append(position)
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
        sites[where_to_append]["idx"] = np.append(sites[where_to_append]["idx"], sites["vacant"]["idx"][i])

    # Delete the ith element from 'vacant'
    sites["vacant"]["positions"] = np.delete(sites["vacant"]["positions"], i_blocked, axis=0)
    sites["vacant"]["translation_vectors"] = np.delete(sites["vacant"]["translation_vectors"], i_blocked, axis=0)  
    sites["vacant"]["idx"] = np.delete(sites["vacant"]["idx"], i_blocked)  
    # For shift_vectors and distances, delete both rows and columns at the indices in i_blocked
    sites["vacant"]["shift_vectors"] = np.delete(sites["vacant"]["shift_vectors"], i_blocked, axis=0)
    sites["vacant"]["shift_vectors"] = np.delete(sites["vacant"]["shift_vectors"], i_blocked, axis=1)
    
    sites["vacant"]["distances"] = np.delete(sites["vacant"]["distances"], i_blocked, axis=0)
    sites["vacant"]["distances"] = np.delete(sites["vacant"]["distances"], i_blocked, axis=1)

def fill_cell_with_closest_hydrogen(sites, blocking_radius):
    n_vacant_sites = len(sites["vacant"]["positions"])
    t_now = time.time()
    iteration = 0
    while n_vacant_sites != 0:
        if iteration == 0:    
            ith_atom = np.random.randint(0, n_vacant_sites)
        else:
            # Calculate the distances between occupied and vacant sites
            distances = np.linalg.norm(sites["occupied"]["positions"][:, np.newaxis] - sites["vacant"]["positions"], axis=2).min(axis=0)
            # Filter out distances that are less than the blocking radius by setting them to infinity (or another large value)
            valid_distances = np.where(distances <= blocking_radius, distances, np.inf)
            # Get the index of the smallest valid distance
            ith_atom = np.argmin(valid_distances)
        #if iteration % 10 == 0:
        #    print(f"Time for {iteration} / {n_vacant_sites} possible iterations: ", time.time() - t_now)
        iteration += 1
        # Get the indices where the distance is less than the threshold
        indices = np.where(sites["vacant"]["distances"][ith_atom] < blocking_radius)[0]
        delete_and_append_for_matrix(sites, np.concatenate(([ith_atom], indices)), ith_atom)
        n_vacant_sites -= len(indices) + 1
    return sites
