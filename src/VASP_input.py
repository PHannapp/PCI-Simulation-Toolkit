from bson.objectid import ObjectId  # Required to work with _id
import os
import numpy as np
import gridfs
from src.CONSTANTS import MAGNETIZATION_DICT
import shutil
from datetime import datetime, timezone  # Import timezone for UTC handling
from MongoDB.connect import collection_calculation, collection_system, db

# Save the current working directory
original_directory = os.getcwd()

directory_mcsqs = os.path.join(original_directory, "RUNS")
directory_vasp_results_base = os.path.join(original_directory, "VASPResults")
directory_vasp_results =  os.path.join(directory_vasp_results_base, "copy")
directory_script = os.path.join(original_directory, "VASPInput")
directory_vasp_input = os.path.join(original_directory, "VASPInput")
POTCAR_dir = os.path.join(original_directory,  "POTCARS")
target_directory =  os.path.join(directory_vasp_results_base, "old")

def vasp_input():
    tasks = collection_calculation.find(
        {
            "$or": [
                {'POSCAR': {'$exists': True, '$ne': None}, 'TO_RUN':False}
            ]
            
        }
    )
    EDIFF_BASE = [1e-5]
    EDIFFG = [1e-4]
    LREAL = ["Auto"]
    SAVE = ['.TRUE.']
    ISMEAR = ['']
    ISIF = [7]
    for task in tasks:
        if not task.get('system_id'):
            task_system = {'calculation_information': {'cutoff':[450], 'k_point_density': [0.22]}}
        else:   
            task_system = collection_system.find_one({"_id": task['system_id']})

        # Create folder name
        task_id_str = ""  
        for elements in task['elements']:
            for element in elements:
                task_id_str += element
            task_id_str += "-"
        task_id_str += f'{str(task["_id"])}'  # Convert ObjectId to string
        task_dir = os.path.join(directory_vasp_input, task_id_str)  # Create the path for the folder

        # Check if the directory_VASP already exists; if not, create it
        if not os.path.exists(task_dir):
            os.makedirs(task_dir)
            print(f"Created directory: {task_dir}")
        poscar_array = task["POSCAR"].splitlines()
        elements = poscar_array[5].split()
        number = poscar_array[6].split()
        number = [int(n) for n in number]

        with open(os.path.join(directory_script,"script.sh"), 'r') as script_file:
            # Read the content of the individual POTCAR file and write it to the combined POTCAR file
            script_string = script_file.read()
        # Create the POTCAR file by concatenating all element-specific POTCAR files
        with open(os.path.join(task_dir, 'POTCAR'), 'w') as POTCARfile:
            potcar_string = ""
            for element in elements:
                potcar_path = os.path.join(POTCAR_dir, element.upper())  # Find the path to the POTCAR file for each element
                # Ensure that the POTCAR file exists before reading
                if os.path.exists(potcar_path):
                    with open(potcar_path, 'r') as potcar_element_file:
                        # Read the content of the individual POTCAR file and write it to the combined POTCAR file
                        potcar_string_tmp = potcar_element_file.read()
                        potcar_string += potcar_string_tmp
                        POTCARfile.write(potcar_string_tmp)
                else:
                    print(f"POTCAR for element {element.upper()} not found in {POTCAR_dir}")

        # Create the POSCAR file 
        with open(os.path.join(task_dir, 'POSCAR'), 'w') as POSCARfile:
                POSCARfile.write(task["POSCAR"])
        with open(os.path.join(task_dir, 'script.sh'), 'w') as SCRIPTfile:
                SCRIPTfile.write(script_string)
        # Extract the 2nd to 4th lines from the POSCAR file and calculate the lattice parameters
        lattice_vectors = [np.array(list(map(float, poscar_array[i].split()))) for i in range(2, 5)]
        # Calculate the magnitudes of the lattice vectors (a, b, c)
        lattice_params = [np.linalg.norm(vec) for vec in lattice_vectors]
        k_point_string2 = None
        k_point_density = task_system["calculation_information"]["k_point_density"][0]
        # Get the k-point density from the task dictionary
        # Calculate the k-point sampling for each lattice parameter
        kpoints = [int(np.ceil(2 * np.pi / (lp * k_point_density))) for lp in lattice_params]
        # Create the KPOINTS file 
        with open(os.path.join(task_dir, 'KPOINTS'), 'w') as KPOINTSfile:
            k_point_string = f"""Automatic
0
Gamma
{kpoints[0]} {kpoints[1]} {kpoints[2]}
    """
            KPOINTSfile.write(k_point_string)
        system = ""
        mag = ""
        for i in range(len(elements)):
            if elements[i].upper() == "VA":
                continue
            system += elements[i] + str(number[i])
            mag += f"{str(number[i])}*{MAGNETIZATION_DICT[elements[i].upper()]} "
        EDIFF = EDIFF_BASE[0] #* sum(number)
        ISPIN = ""
        ISPIN = f"#ISPIN=2\n#MAGMOM={mag}"
        # Write "HI" into the INCAR file
        with open(os.path.join(task_dir, 'INCAR'), 'w') as incar_file:
            incar_string = f"""SYSTEM =  {system}
Accuracy and Algorithm:
PREC= Normal 
IALGO=38   # RMM-DIIS
NPAR=2
ENCUT={task_system['calculation_information']['cutoff'][0]} 
EDIFF = {EDIFF}
LWAVE={SAVE[0]}
LCHARG={SAVE[0]}

 Ionic Relaxation:
NSW= 200
EDIFFG = {EDIFFG[0]}
IBRION=2   # 1 RMM-DIIS, 2 CG for ionic relaxation
ISIF={ISIF[0]}
LREAL={LREAL[0]}
{ISPIN}
Smearing and DOS related values:
{ISMEAR[0]}
SIGMA = 0.2
"""
            incar_file.write(incar_string)

        # Build the $set dictionary dynamically
        update_set = {
            f"VASP.INCAR": incar_string,
            f"VASP.KPOINTS": k_point_string,
            f"VASP.POSCAR": task["POSCAR"],
            f"VASP.POTCAR": potcar_string,
            f"VASP.SCRIPT": script_string,
            "folder_str": task_id_str,
            "TO_RUN": True
        }

        # Perform the update operation
        collection_calculation.update_one(
            {"_id": task["_id"]},
            {
                "$set": update_set
            }
        )

fs = gridfs.GridFS(db)  # Initialize GridFS

def process_output():
    # Iterate over each folder in the directory
    for folder in os.listdir(directory_vasp_results):
        folder_path = os.path.join(directory_vasp_results, folder)
        print(f'Opening folder {folder}')
        if folder =="old":
            continue
        # Ensure it's a folder
        if os.path.isdir(folder_path):
            try:
                # Extract the 'id' and 'source' from the folder name
                id_str = folder.split('-')[-1]
                task = collection_calculation.find_one({"_id": ObjectId(id_str)})
                if not task:
                    print(f"Task with id {id_str} not found in the database.")
                    continue
                collection_calculation.update_one(
                    {"_id": task["_id"]},
                    {"$set": {f"timestamp_analyzed": datetime.now(timezone.utc)}}
                )
                # Update the database with the file content under 'VASP.{source}.{filename}'
                #collection_calculation.update_one(
                #    {"_id": task["_id"]},
                #    {"$unset": {f"VASP.{source}": ""}}
                #)
                energy = None
                forces = None
                file_name = None
                total_forces = None
                # Iterate over each file in the folder
                for file_name in os.listdir(folder_path):
                    ALREADY_STORED = False
                    if file_name.endswith(".h5") or file_name.endswith(".sh"):
                        print(f"Skipping binary file: {file_name}")
                        continue
                    file_path = os.path.join(folder_path, file_name)
                    
                    # Skip files that are too large unless it's OUTCAR
                    file_size = os.path.getsize(file_path)  # Get file size in bytes
                    max_field_size = 16777216  # MongoDB's default BSON document limit (16MB)
                    
                    # Read file content if it doesn't exceed the size limit, unless it's OUTCAR
                    if file_size > max_field_size / 10 and 'OUTCAR' not in file_name:
                        print(f"Skipping {file_name} as it's too large.")
                        continue

                    base_file_name = os.path.basename(file_name).upper()
                    
                    if file_size > max_field_size / 10 and 'OUTCAR' in file_name:
                        print(f"Chunking large OUTCAR file into GridFS: {file_name}")
                        
                        # Store the OUTCAR file in GridFS
                        with open(file_path, 'rb') as file:  # Open in binary mode
                            # Upload the file to GridFS and get the file ID
                            grid_outcar_id = fs.put(file, filename=file_name)
                        
                        # Update MongoDB to store the GridFS ID reference
                        collection_calculation.update_one(
                            {"_id": task["_id"]},
                            {"$set": {f"VASP.{base_file_name}": grid_outcar_id}}
                        )
                        print(f"OUTCAR file stored in GridFS with ID: {grid_outcar_id}")
                        ALREADY_STORED = True
                        # Read file content
                    if file_name == "OUTCAR":
                        with open(file_path, 'r') as file:
                            file_content = file.read()
                        energy_outcar, forces_outcar, total_forces = extract_energy_force_outcar(file_content)
                        energy = energy_outcar if energy_outcar else energy
                        forces = forces_outcar if forces_outcar else forces
                        converged = check_convergence_outcar(file_content)

                        if total_forces :
                            if total_forces  <= 0.05:
                                converged = True
                            else:
                                converged = False
                        else:
                            converged = False

                        # Update MongoDB with the total forces from OUTCAR
                        collection_calculation.update_one(
                            {"_id": task["_id"]},
                            {"$set": {f"VASP.total_forces": total_forces,
                                      f"VASP.energy": energy,
                                      f"VASP.forces": forces,
                                      f"VASP.SUCCESS": converged}}
                        )
                        continue
                    if not ALREADY_STORED:
                        # Read file content
                        with open(file_path, 'r') as file:
                            file_content = file.read()

                        # Determine file type (e.g., POSCAR, CONTCAR, OUTCAR, vasprun.xml)

                        ## Check for OUTCAR or vasprun.xml
                        #if base_file_name == "OUTCAR":
                        #    energy_outcar, forces_outcar = extract_energy_force_outcar(file_content)
                        #    energy = energy_outcar if energy_outcar else energy
                        #    forces = forces_outcar if forces_outcar else forces
                        #    #converged = check_convergence_outcar(file_content)
                        #    total_forces = extract_total_forces_oszicar(file_content)
#   
                        #    # Update MongoDB with the total forces from OSZICAR
                        #    collection_calculation.update_one(
                        #        {"_id": task["_id"]},
                        #        {"$set": {f"VASP.{source}.total_forces": total_forces}}
                        #    )

                        #if base_file_name == "VASPRUN.XML":
                        #    energy_vasp, forces_vasp = extract_energy_force_vasprun(file_content)
                        #    energy = energy_vasp if energy_vasp else energy
                        #    forces = forces_vasp if forces_vasp else forces
                        #    #converged = check_convergence_vasprun(file_content)
                        #
                        # Update the database with the file content under 'VASP.{source}.{filename}'
                        collection_calculation.update_one(
                            {"_id": task["_id"]},
                            {"$set": {f"VASP.{base_file_name}": file_content}}
                        )
                                
                # Prepare destination path
                destination = os.path.join(target_directory, folder)

                # If the destination folder exists, remove it
                if os.path.exists(destination):
                    shutil.rmtree(destination)
                shutil.move(folder_path, target_directory)
                print(f"Moved folder {folder} to {target_directory}")
            except ValueError as e:
                # Prepare destination path
                # Move folder to target directory even if there is an error
                if folder_path != target_directory:
                    print(f"Error processing folder {folder}/{file_name}: {e}")
                    destination = os.path.join(target_directory, folder)

                    # If the destination folder exists, remove it
                    if os.path.exists(destination):
                        shutil.rmtree(destination)
                    shutil.move(folder_path, target_directory)
                    print(f"Moved folder {folder} to {target_directory} after error")


def check_convergence_outcar(content):
    """
    Check if OUTCAR indicates convergence.
    """
    # Check for a string indicating convergence in OUTCAR
    if "reached required accuracy" in content.lower():
        return True
    return False

def extract_total_forces_oszicar(content):
    """
    Extract the total forces from the OSZICAR file.
    """
    total_forces = None
    # Split the content into lines and reverse it to get the last occurrence of 'FORCES:'
    for line in reversed(content.splitlines()):
        if "FORCES:" in line:
            total_forces = float(line.split()[-1])  # Extract the last value in the line as float
            break
    return total_forces

def check_convergence_vasprun(content):
    """
    Check if vasprun.xml indicates convergence.
    """
    # Look for <calculation> and check if electronic and ionic steps are converged
    if "<calculation>" in content and "<converged>True</converged>" in content:
        return True
    return False

def extract_energy_force_outcar(content):
    """
    Extract energy and forces from the last appearance of the TOTAL-FORCE section in OUTCAR file.
    """
    energy = None
    forces = []
    max_force_component = None

    # Extract the total energy (TOTEN)
    for line in content.splitlines():
        if "free  energy   TOTEN" in line:
            energy = float(line.split()[-2])  # Extract the energy value

    # Reverse iterate to find the last TOTAL-FORCE section
    lines = content.splitlines()
    force_section_found = False
    for i in range(len(lines) - 1, -1, -1):
        if "TOTAL-FORCE" in lines[i]:
            force_section_found = True
            force_start_index = i + 1  # Start extracting forces after this line
            break

    # If TOTAL-FORCE section is found, extract the forces
    if force_section_found:
        for line in lines[force_start_index+1:]:
            if "---" in line:  # End of the TOTAL-FORCE section
                break
            forces_tmp = [float(val) for val in line.split()[3:6]]
            forces.append(forces_tmp)  # Extract force vectors
            max_component_in_line = max(abs(comp) for comp in forces_tmp)
            if max_force_component is None or max_component_in_line > max_force_component:
                max_force_component = max_component_in_line

    return energy, forces, max_force_component

import xml.etree.ElementTree as ET

def extract_energy_force_vasprun(content):
    """
    Extract energy and forces from vasprun.xml file.
    """
    energy = None
    forces = []

    # Check if the file ends with the correct closing tag
    if not content.strip().endswith("</modeling>"):
        print("Error: vasprun.xml seems to be incomplete.")
        return energy, forces

    try:
        # Parse the XML content
        root = ET.fromstring(content)

        # Find the forces
        for varray in root.findall(".//varray[@name='forces']"):
            for v in varray.findall("v"):
                forces.append([float(f) for f in v.text.split()])

        # Extract the final energy from the last calculation step
        for scstep in root.findall(".//calculation/energy/i[@name='e_fr_energy']"):
            energy = float(scstep.text)  # The free energy from the final calculation step

    except ET.ParseError as e:
        print(f"Error parsing vasprun.xml: {e}")

    return energy, forces