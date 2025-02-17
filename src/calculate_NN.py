import numpy as np

def seperate_vac_atoms(atoms):
    vac_atoms = []
    for atom in atoms:
        if "Vac" in atom.label:
            vac_atoms.append(atom)
    return vac_atoms

def delete_and_append(sites, i_blocked, i_occupied):
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
    sites["vacant"]["idx"] = np.delete(sites["vacant"]["idx"], i_blocked)

def calculate_nearest_neighbors(lattice_vectors, vac_atom, vac_atoms, maximum_distance,verbose=False):
        atom_position = np.array([float(pos) for pos in vac_atom.position])
        if verbose:
            print("Atom: \n", vac_atom.label)
            print("atom_position: \n", atom_position)
        # Check how many translation vectors in each lattice direction, to reach 5 A
        translation_vector_bounds = np.zeros((3, 2))
        for axis, lattice_vector in enumerate(lattice_vectors):
            for direction in range(2):
                n = 0
                if direction == 0: # negative direction
                    while (atom_position[axis] + n) * np.linalg.norm(lattice_vector) < maximum_distance:
                        n += 1
                elif direction == 1: # positive direction
                    while ((1 - atom_position[axis]) + n) * np.linalg.norm(lattice_vector) < maximum_distance:
                        n += 1
                translation_vector_bounds[axis, direction] = n
        if verbose:
            print("translation_vector_bounds: \n", translation_vector_bounds)
        # Adjust for both negative and positive translations
        x_translation = np.arange(-translation_vector_bounds[0][0], translation_vector_bounds[0][1] + 1) 
        y_translation = np.arange(-translation_vector_bounds[1][0], translation_vector_bounds[1][1] + 1) 
        z_translation = np.arange(-translation_vector_bounds[2][0], translation_vector_bounds[2][1] + 1)   
        # Use meshgrid to get all combinations of x, y, z translations
        x, y, z = np.meshgrid(x_translation, y_translation, z_translation, indexing='ij')

        # Combine the meshgrid into a single array of translation vectors
        all_translation_vectors = np.vstack([x.ravel(), y.ravel(), z.ravel()]).T
        amount_of_unit_cells = len(all_translation_vectors)
        # create new atom array, where a copy of the original and translated atoms are inside
        all_distances = np.zeros((len(vac_atoms) * amount_of_unit_cells)) 
        # Fill the new atom array
        for ith_unit_cell, translation_vector in enumerate(all_translation_vectors):
            add = len(vac_atoms) * ith_unit_cell
            for nth_atom_in_original_unit_cell in range(len(vac_atoms)):
                # position of translated atom
                other_atom_position = np.array([float(pos) for pos in vac_atoms[nth_atom_in_original_unit_cell].position]) + translation_vector

                # calculate distance of mth_atom to atom_position
                dij = np.linalg.norm(np.dot(lattice_vectors, other_atom_position - atom_position))

                # append to all_distances
                mth_atom = nth_atom_in_original_unit_cell + add
                all_distances[mth_atom] = dij

        n_nearest_neighbors = np.sum(all_distances < maximum_distance) - 1 # -1 because the original atoms has distance 0 to itself
        if verbose:
            print("Number of nearest neigbors: ", n_nearest_neighbors)
        return max(n_nearest_neighbors, 0)
