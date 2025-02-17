import numpy as np

#calcualte error:
def calculate_supercell_volume(lattice_vectors, supercell_size):
    """
    Calculate the volume of the supercell based on the lattice vectors and supercell size.
    """
    # Calculate the volume of the unit cell (scalar triple product)
    unit_cell_volume = np.abs(np.dot(lattice_vectors[0], np.cross(lattice_vectors[1], lattice_vectors[2])))
    
    # Supercell volume is scaled by supercell_size^3
    print(unit_cell_volume)
    supercell_volume = unit_cell_volume * (supercell_size ** 3)
    return supercell_volume

def calculate_surface_area(lattice_vectors, supercell_size):
    """
    Calculate the surface area of the supercell based on the lattice vectors and supercell size.
    """
    # Get the magnitudes of the lattice vectors
    a1 = np.linalg.norm(lattice_vectors[0])
    a2 = np.linalg.norm(lattice_vectors[1])
    a3 = np.linalg.norm(lattice_vectors[2])

    # Each face has two vectors that form it. For example, for one face:
    # Surface area of face formed by a1 and a2 is |a1 x a2|
    surface_area_face1 = a1 * a2 * supercell_size**2
    surface_area_face2 = a2 * a3 * supercell_size**2
    surface_area_face3 = a1 * a3 * supercell_size**2

    # Total surface area (sum of all 6 faces)
    total_surface_area = 2 * (surface_area_face1 + surface_area_face2 + surface_area_face3)
    return total_surface_area

def calculate_blocked_volume(blocking_radius, surface_area):
    """
    Calculate the blocked volume near the surface based on the blocking radius and surface area.
    """
    # Blocked volume is approximately the surface area * blocking_radius (for thin blocked regions near the surface)
    blocked_volume = surface_area * blocking_radius
    return blocked_volume

def estimate_error(supercell_volume, blocked_volume):
    """
    Estimate the error based on the blocked volume outside the supercell and the total volume.
    """
    # Error is the ratio of blocked volume to the total volume of the supercell
    error_estimate = blocked_volume / supercell_volume
    return error_estimate

def estimate_error_sweep(supercell_size, lattice_vectors, blocking_radius):
    # Step 1: Calculate the supercell volume
    supercell_volume = calculate_supercell_volume(lattice_vectors, supercell_size)
    print(f"Supercell Size: {supercell_size}")
    # Step 2: Calculate the surface area of the supercell
    surface_area = calculate_surface_area(lattice_vectors, supercell_size)
    # Step 3: Calculate the blocked volume due to the blocking radius
    blocked_volume = calculate_blocked_volume(blocking_radius, surface_area)
    # Step 4: Estimate the error
    error = estimate_error(supercell_volume, blocked_volume)
    print(f"Estimated Error: {error * 100:.2f}%")