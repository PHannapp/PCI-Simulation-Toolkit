import os
from ase.io.vasp import read_vasp_out
import matplotlib.pyplot as plt
import imageio.v2 as imageio
from ase.io import write
import tempfile
import numpy as np
from PIL import Image

# Input and output paths
outcar_path = "/home/pet69516/Dokumente/DFTforCALPHAD/VASPResults/old/La-AlNi-H-Va-671dba1ac2ff5af1b377549c_ISIF7/OUTCAR_SAVE1"
output_gif_path = "/home/pet69516/Dokumente/DFTforCALPHAD/CALPHAD_CALC/OUTCAR_gif/"

# Parameters for customization
rotation_angle = '80z, 00y, 100x'
atom_radii_dict = {'H': 0.4, 'Ti': 1, 'Mn': 0.8}

# Extracting atomic structures from OUTCAR file
atoms_list = list(read_vasp_out(outcar_path, index=slice(None)))

# Create a temporary directory to store frames
with tempfile.TemporaryDirectory() as temp_dir:
    frames = []
    max_size = None

    # Loop through all steps and generate an image for each step
    for i, atoms in enumerate(atoms_list):
        radii = [atom_radii_dict.get(atom.symbol, 1.0) for atom in atoms]

        frame_path = os.path.join(temp_dir, f'frame_{i}.png')
        
        # Use ASE to visualize and save the current atomic configuration
        fig, ax = plt.subplots(dpi=300)  # Increase the DPI for higher resolution
        ax.axis('off')

        write(frame_path, atoms, show_unit_cell=2, rotation=rotation_angle, radii=radii)
        plt.close()

        # Open the image and resize all frames to the same size
        img = Image.open(frame_path)
        if max_size is None:
            max_size = img.size
        else:
            img = img.resize(max_size, Image.Resampling.LANCZOS)

        frames.append(np.array(img))

    name = f"{outcar_path.split('/')[-2]}.gif"
    # Combine all frames into a GIF
    imageio.mimsave(os.path.join(output_gif_path, name), frames, fps=5)

print(f"GIF saved successfully at: {output_gif_path}/{name}")
