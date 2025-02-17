from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Load the POSCAR file
poscar = Poscar.from_file("/home/pet69516/Dokumente/DFTforCALPHAD/POSCARS/TEST_POSCAR_LaAlVaHVa_after")

# Create structure object
structure = poscar.structure

# Analyze symmetry
analyzer = SpacegroupAnalyzer(structure)

# Get space group and symmetry information
space_group = analyzer.get_space_group_symbol()
crystal_system = analyzer.get_crystal_system()
symmetry_dataset = analyzer.get_symmetry_dataset()

# Output the results
print(f"Space group: {space_group}")
print(f"Crystal system: {crystal_system}")
#print("Symmetry dataset:")
#print(symmetry_dataset)