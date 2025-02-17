from pymongo import MongoClient
from bson.objectid import ObjectId  # Required to work with _id
import numpy as np
from CONSTANTS import metallic_radii_dict

# Connect to MongoDB
client = MongoClient('mongodb://localhost:27017/')
db = client.DFTforCALPHAD
collection_system = db.System
collection_calculation = db.calculation

def write_poscar(endmember, system_id):
    task_system = collection_system.find_one({"_id":ObjectId(system_id)})
    # check if exception
    poscar_string = ""
    if not task_system["exceptions"]:
        structure = task_system["structure"]

    else:
        NOBREAK = False
        for exception in task_system["exceptions"]:
            if NOBREAK: 
                break
            NOBREAK = True
            structure = exception["structure"]
            for element in range(len(endmember)):
                if exception["elements"][element] != [] and endmember[element] != exception["elements"][element]:
                    NOBREAK = False
                    structure = task_system["structure"]
                    break
    number_of_atoms = [len(site) for site in structure["atomic_positions"]]
    # Write the file line by line
    # 1. line: Comment line 
    line  = ""
    for element in range(len(number_of_atoms)):
        line += f"{endmember[element][0]}{number_of_atoms[element] }"
    poscar_string += line + "\n"
    # 2. Scaling factor:
    poscar_string += "1.0\n"
    multiplicator = 0
    devisor = 0
    for site in range(len(task_system["reference_solids"])):
        multiplicator += metallic_radii_dict[task_system["reference_solids"][site].lower()] / metallic_radii_dict[endmember[site][0].lower()] * task_system["multiplicity"][site]
        devisor += task_system["multiplicity"][site]
    
    # 3. line: lattice vectors
    for vectors in structure["lattice_parameters"]:
        for vector in vectors:
            poscar_string += str(float(vector) / multiplicator * devisor) + ' '
        poscar_string += "\n"
    # 4. line: type and number of atoms
    for element in range(len(number_of_atoms)):
        poscar_string += f"{endmember[element][0]} "
    poscar_string += "\n"
    for element in range(len(number_of_atoms)):
        poscar_string += f"{number_of_atoms[element]} "
    poscar_string += "\ndirect\n"
    # 5. line: Atomic positions
    for atom_type in structure["atomic_positions"]:
        for atomic_positions in atom_type:
            for position in atomic_positions:
                poscar_string += str(position) + ' '
            poscar_string += "\n"
    return poscar_string

#def write_mcsqs_inputs():
#    task_systems = collection_system.find({})
#    for system in task_systems:
#        for site in range(len(system["elements"])):
#            for proportion in [0.125,0.25,0.375,0.5]:


def write_rndstr():
    tasks = collection_calculation.find(
        {
            'type': 'mcsqs',
            'RNDSTR': {'$exists': False}
        }
    )

    for task in tasks:
        task_system = collection_system.find_one(
            {
                '_id': task['system_id']
            }) 
        # Initialize an array of dictionaries for each sublattice
        endmembers = [[], []]
        for i, sublattice in enumerate(task['elements']):
            if len(sublattice) == 2:
                for sl in range(2):
                    endmembers[sl].append([sublattice[sl]])
            else:
                for sl in range(2):
                    endmembers[sl].append(sublattice)
        lattice_parameters = [[], []]
        break_occurred = False  # Set the flag if a break occurs
        for sl, endmember in enumerate(endmembers):
            task_em = collection_calculation.find_one(
                {
                    'elements': endmember,
                    'system_id': task['system_id'],
                    'RNDSTR': {'$exists': False}
                }
            )
            if not (task_em.get('VASP') and task_em['VASP'].get('CONTCAR')):
                break_occurred = True  # Set the flag if a break occurs
                break
            else:
                vectors_tmp = []
                contcar_array = task_em['VASP']['CONTCAR'].splitlines()
                for i in range(3):
                    vectors = contcar_array[2 + i].split()
                    vectors_tmp.append([float(vec) for vec in vectors])
                    lattice_parameters[sl].append((vectors_tmp[i][0]**2 + vectors_tmp[i][1]**2 + vectors_tmp[i][2]**2)**0.5)
        if break_occurred:
            continue
        ratio = [round(lattice_parameters[0][i] / lattice_parameters[1][i], 3) for i in range(3)]
        mixing_species = task["elements"]
        # check if exception
        rndstr_string = ""
        if not task_system["exceptions"]:
            structure = task_system["structure"]
            unit_cell_size = task_system["unit_cell_size"]
        else:
            NOBREAK = False
            for exception in task_system["exceptions"]:
                if NOBREAK: 
                    break
                NOBREAK = True
                structure = exception["structure"]
                unit_cell_size = exception["unit_cell_size"]
                for element in range(len(mixing_species)):
                    if exception["elements"][element] != [] and mixing_species[element] != exception["elements"][element]:
                        NOBREAK = False
                        structure = task_system["structure"]
                        unit_cell_size = task_system["unit_cell_size"]
                        break
        number_of_atoms = [len(site) for site in structure["atomic_positions"]]
        # Write the file line by line   
        # 1. line: lattice vectors
        proportion = task['mixing_proportion']
        for sl, vectors in enumerate(vectors_tmp):
            for vector in vectors:
                rndstr_string += str(vector * ((ratio[sl]-1)*proportion + 1)) + ' '
            rndstr_string += "\n"
        # 2. line: unit cell repetition
        rndstr_string += "1 0 0\n0 1 0\n0 0 1\n"
        # 4. line: type and number of atoms
        # 5. line: Atomic positions
        for site, atom_type in enumerate(structure["atomic_positions"]):
            if len(mixing_species[site]) == 2:
                line_ending = f"{mixing_species[site][0]}={round(proportion, 3)}, {mixing_species[site][1]}={round(1 - proportion, 3)}\n"
            else:
                line_ending = f"{mixing_species[site][0]}\n"
            for atomic_positions in atom_type:
                for position in atomic_positions:
                    rndstr_string += str(position) + ' '
                rndstr_string += line_ending
        sqscell_string = f"""1

{unit_cell_size[0]} 0 0
0 {unit_cell_size[1]} 0
0 0 {unit_cell_size[2]}"""
        collection_calculation.update_one({'_id': task['_id']}, {'$set': {'RNDSTR': rndstr_string, 'SQSCELL': sqscell_string, 'TO_RUN': False}})
        print(f"Created RNDSTR for {task['elements']}")
        
        #return rndstr_string

def write_mcsqs(task):
    task_system = collection_system.find_one({"_id": ObjectId(task["system_id"])})
    # Calculate the length of the range of nearest neighbors, to calcualte the correlation function
    lattice_vectors = task_system["structure"]["lattice_parameters"]
    shift_range = range(-1, 2)  # Shifts from -1 to 1
    # Create a grid of neighboring points to calculate the distances
    final_neighbors = []
    for n1 in shift_range:
        for n2 in shift_range:
            for n3 in shift_range:
                if (n1, n2, n3) != (0, 0, 0):  # Ignore the origin
                    neighbor_vector = n1 * np.array(lattice_vectors[0]) + n2 * np.array(lattice_vectors[1]) + n3 * np.array(lattice_vectors[2])
                    distance = np.linalg.norm(neighbor_vector)
                    final_neighbors.append(distance)

    # Sort the distances to find the third nearest neighbor
    final_neighbors.sort()
    # Display the first few neighbors for context
    final_neighbor_distance = final_neighbors[5]
    # Calculate the size of supercell by looking at total numbers of atoms and number of mixing atoms
    # Minimum unit cell size:
    #for sublattice, site in enumerate(task["elements"]):
    #    if len(site) == 2:
    #        break
#
#
    #if not task_system["exceptions"]:
    #    unit_cell_size = task_system["unit_cell_size"]
    #else:
    #    NOBREAK = False
    #    for exception in task_system["exceptions"]:
    #        if NOBREAK: 
    #            break
    #        NOBREAK = True
    #        unit_cell_size = exception["unit_cell_size"]
    #        for element in range(len(task["elements"])):
    #            if exception["elements"][element] != [] and task["elements"][element] != exception["elements"][element]:
    #                NOBREAK = False
    #                unit_cell_size = task_system["unit_cell_size"]
    #                break
#
#
    #minimum_unit_cell_size = unit_cell_size[sublattice]
    #number_of_atoms = len(task["RNDSTR"].splitlines()) - 6
    #lines_with_char = [line for line in task["RNDSTR"].splitlines() if "=" in line]
    #number_of_mixing_atoms = len(lines_with_char)
    #while True:
    #    if number_of_mixing_atoms * minimum_unit_cell_size * task["mixing_proportion"] % 1 == 0:
    #        size_of_supercell = number_of_atoms * minimum_unit_cell_size
    #        break
    #    minimum_unit_cell_size += 1    
#
    mcsqs_string = f"""
import subprocess
import time
import psutil

def kill_process_and_children(pid):
    try:
        parent = psutil.Process(pid)
        children = parent.children(recursive=True)  # Get child processes
        for child in children:
            child.kill()  # Kill child processes
        parent.kill()  # Kill the parent process
        parent.wait()  # Wait for the parent to be fully terminated
    except psutil.NoSuchProcess:
        print("Process ", pid, " does not exist")

# Open the log file in append mode so output can be written continuously
print("Starting MCSQS calculation")
final_neighbor_distance = {final_neighbor_distance}
cmd_string = "/home/pet69516/ATAT/bin/mcsqs -2=" + str(final_neighbor_distance) + " -3=" + str(final_neighbor_distance)
with open("stdout.txt", "a") as stdout_file, open("stderr.txt", "a") as stderr_file:
    # Start the process and redirect stdout and stderr to the file
    process = subprocess.Popen(cmd_string,
                               shell=True, stdout=stdout_file, stderr=stderr_file, bufsize=1, universal_newlines=True)
    
    while process.poll() is None:  # While the process is running
        stdout_file.flush()  # Force the file to be updated
        stderr_file.flush()  # Force the file to be updated
        # Optionally, add a short sleep to avoid tight looping
        time.sleep(10)

    # Final flush to ensure all remaining output is written
    stdout_file.flush()  # Force the file to be updated
    stderr_file.flush()  # Force the file to be updated
    process.wait()
    kill_process_and_children(process.pid)
print("First run finished")
# Reopen the stdout.txt file to count lines
with open("stdout.txt", "r") as stdout_file:
    stdout_content = stdout_file.read()  # Read the content of the file
    stdout_lines = stdout_content.splitlines()  # Split the content into lines

while len(stdout_lines) < 25:
    final_neighbor_distance += 1
    print("No result of first step (Number of lines: ", len(stdout_lines), "). New neighbour distance:", final_neighbor_distance)
    cmd_string = "/home/pet69516/ATAT/bin/mcsqs -2=" + str(final_neighbor_distance) + " -3=" + str(final_neighbor_distance)
    with open("stdout.txt", "w") as stdout_file, open("stderr.txt", "a") as stderr_file:
        # Start the process and redirect stdout and stderr to the file
        process = subprocess.Popen(cmd_string,
                                   shell=True, stdout=stdout_file, stderr=stderr_file, bufsize=1, universal_newlines=True)

        while process.poll() is None:  # While the process is running
            stdout_file.flush()  # Force the file to be updated
            stderr_file.flush()  # Force the file to be updated
            # Optionally, add a short sleep to avoid tight looping
            time.sleep(10)

        # Final flush to ensure all remaining output is written
        stdout_file.flush()  # Force the file to be updated
        stderr_file.flush()  # Force the file to be updated 
    # Reopen the stdout.txt file to count lines
    with open("stdout.txt", "r") as stdout_file:
        stdout_content = stdout_file.read()  # Read the content of the file
        stdout_lines = stdout_content.splitlines()  # Split the content into lines
    process.wait()    
    kill_process_and_children(process.pid)

t_start = time.time()
with open("stdout_main.txt", "a") as stdout_file, open("stderr_main.txt", "w") as stderr_file:
    # Start the process and redirect stdout and stderr to the file
    process = subprocess.Popen(f"/home/pet69516/ATAT/bin/mcsqs -rc",
                               shell=True, stdout=stdout_file, stderr=stderr_file, bufsize=1, universal_newlines=True)
    
    while process.poll() is None:  # While the process is running
        if (time.time() - t_start) >= 60 * 60 * 0.2: # five hour
            kill_process_and_children(process.pid)
            time.sleep(5)
            if process.poll() is None:
                kill_process_and_children(process.pid)
            break
        stdout_file.flush()  # Force the file to be updated
        stderr_file.flush()  # Force the file to be updated
        # Optionally, add a short sleep to avoid tight looping
        time.sleep(10)

    # Final flush to ensure all remaining output is written
    stdout_file.flush()  # Force the file to be updated
    stderr_file.flush()  # Force the file to be updated

process = subprocess.Popen(f"sqs2poscar bestsqs.out", shell=True)    
process.wait()   
kill_process_and_children(process.pid)
time.sleep(5)
    """
    return mcsqs_string

def clean_mcsqs_poscar(uncleaned):

    # Split the content into individual lines
    uncleaned_lines = uncleaned.splitlines()
    
    # Extract the elements (line 5) and the number of atoms (line 6) from the POSCAR
    elements = uncleaned_lines[5].split()
    number_of_atoms = uncleaned_lines[6].split()

    # Search for the element 'Va' (which usually stands for vacancy) in the element list
    for site, element in enumerate(elements):
        if element == 'Va':  # If 'Va' (vacancy) is found, we stop
            break
    else:
        # If 'Va' is not found, there's no vacancy to clean, return the original POSCAR
        number_of_atoms_ints = [int(i) for i in number_of_atoms]
        return False, uncleaned, sum(number_of_atoms_ints)

    number_of_atoms_ints = [int(i) for i in number_of_atoms]
    # Calculate indices for the positions of atoms before and after the removed element
    before = sum(number_of_atoms_ints[:site]) + 8  # Corrected to convert to int and offset by 7 (lines above coordinates)
    after = sum(number_of_atoms_ints[:site + 1]) + 8  # Line number for the start of the remaining coordinates
    # Remove 'Va' from the elements list and its corresponding atom count
    elements.pop(site)
    number_of_atoms.pop(site)

    # Joining the cleaned elements and atom counts back into single strings
    line_five = " ".join(elements)  # Join cleaned elements into a single line
    line_six = " ".join(number_of_atoms)  # Join cleaned atom counts into a single line

    # Reconstruct the cleaned POSCAR file
    cleaned_poscar = uncleaned_lines[:5]  # First five lines (header and scaling)
    cleaned_poscar.extend([line_five, line_six])  # Add the cleaned element and atom count lines
    cleaned_poscar.extend(uncleaned_lines[7:before])  # Add atom positions before the removed element
    cleaned_poscar.extend(uncleaned_lines[after:])  # Add atom positions after the removed element

    # Join the lines into a final POSCAR string, using newlines to format properly
    poscar_file = "\n".join(cleaned_poscar)
    number_of_atoms_ints = [int(i) for i in number_of_atoms]
    return True, poscar_file, sum(number_of_atoms_ints)  # Return the cleaned POSCAR as a string

#write_mcsqs(collection_calculation.find_one({"type": "mcsqs"}))

def turn_TO_RUN_false_ith(i):
    tasks = collection_calculation.find({
        "ith_sqs": i,
        "TO_RUN":True,
        "type": "mcsqs"})
    for task in tasks:
        collection_calculation.update_one({'_id': task['_id']}, {'$set': {'TO_RUN': False}})
        print("Set to false")
        


