import numpy as np

# Translation vectors for moving the metal atom in all directions (+1 in each axis)
translation_vectors = [
    np.array([0, 0, 0]),   # Translate in the +x direction
    np.array([1, 0, 0]),   # Translate in the +x direction
    np.array([0, 1, 0]),   # Translate in the +y direction
    np.array([0, 0, 1]),   # Translate in the +z direction
    np.array([1, 1, 0]),   # Translate in +x and +y directions
    np.array([1, 0, 1]),   # Translate in +x and +z directions
    np.array([0, 1, 1]),   # Translate in +y and +z directions
    np.array([1, 1, 1]),   # Translate in +x, +y, and +z directions
    np.array([-1, 0, 0]),  # Translate in the -x direction
    np.array([0, -1, 0]),  # Translate in the -y direction
    np.array([0, 0, -1]),  # Translate in the -z direction
    np.array([-1, -1, 0]), # Translate in -x and -y directions
    np.array([-1, 0, -1]), # Translate in -x and -z directions
    np.array([0, -1, -1]), # Translate in -y and -z directions
    np.array([-1, -1, -1]) # Translate in -x, -y, and -z directions
]