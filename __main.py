import json
import time
import subprocess
import time
from pymongo import MongoClient, DESCENDING, ASCENDING
import shutil
import os
import traceback
from write import write_mcsqs, write_rndstr, turn_TO_RUN_false_ith
from read import safe_read_file, safe_read_lines
from create_database_entries_from_system import create_reference_elements, create_endmembers, create_mcsqs
from clean_poscar_from_Va import clean_poscar_from_Va
from VASP_input import vasp_input, process_output, process_output_mcsqs, vasp_input_advanced, directory_vasp_input, directory_vasp_results
from CHGNet_helpers import save_trajectory, CHGNet_calculation
from analyze import print_energies, print_energies_from_OSZICAR, plot_endmember_energies, plot_mixing_energies, show_strcuture, fit_mixing_energies_to_interaction_parameters
from calculate import calculate_energies

# Save the current working directory
original_directory = os.getcwd()
vasp_input()
# modus # TODO -> multiple with 4 or 8 cores simulatneously ... RUNNING
# TODO handle exceptions (if maxstep over)
# TODO magnetic elements:   nspin           = 2     ; For spin-polarized calculation
#  starting_magnetization(1) = 0.0   ; Assuming La is non-magnetic
#  starting_magnetization(2) = 0.05  ; Example value for Ni, adjust based on your system

#system

#execute(testfile, calculation)
#convergence(testfile, calculation, 'Ce')
def read_config():
    default_values = {'THRESHOLD': 0, 'FORCE': 1000}
    try:
        with open('./config.json', 'r') as file:
            config = json.load(file)
            return config.get('THRESHOLD', default_values['THRESHOLD']), config.get('FORCE', default_values['FORCE'])
    except json.JSONDecodeError as e:
        print(f"Error reading config.json: {e}")
        # Optionally, stop the program
        return default_values['THRESHOLD'], default_values['FORCE']
    except FileNotFoundError:
        print("config.json not found. Using default THRESHOLD value.")
        return default_values['THRESHOLD'], default_values['FORCE']
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        # Optionally, stop the program
        return default_values['THRESHOLD'], default_values['FORCE']
    
# Connect to MongoDB
client = MongoClient('mongodb://localhost:27017/')
db = client.DFTforCALPHAD
collection_system = db.System
collection_calculation = db.calculation

processes = {} # dictionary to store running processes
cores_list = {} # total cores in usage
NOTASK = False
while True:
    try:
        THRESHOLD, FORCE = read_config()
        total_cores_in_use = sum(cores_list.values())

        # Forcefully terminate processes if running_processes exceeds FORCE
        while total_cores_in_use > FORCE:
            pid_to_terminate = next(iter(processes))  # Get one PID to terminate
            processes[pid_to_terminate].terminate()  # Terminate the process
            print(f"Forcefully terminating process {pid_to_terminate}.")
            task = collection_calculation.find_one({'RUNNING_pid': pid_to_terminate})
            shutil.rmtree(task['RUNNING_path'])
            collection_calculation.update_one({'RUNNING_pid': pid_to_terminate}, {'$set': {'TO_RUN': True}, '$unset': {'RUNNING_pid': "", 'RUNNING_path': "", 'INPUT_FILE': ""}})
            del processes[pid_to_terminate]  # Remove from tracking
            del cores_list[pid_to_terminate]
            total_cores_in_use = sum(cores_list.values())
            time.sleep(1)  # Give some time for the process to terminate

        if total_cores_in_use < THRESHOLD and total_cores_in_use < FORCE:
            # Find tasks to run    
            task = list(collection_calculation.aggregate([
                {
                    '$match': {
                        '$and': [
                            {'POSCAR': {'$exists': True, '$ne': None}}, 
                            {'TO_RUN':True}, 
                            {'cores': {'$lt': min(THRESHOLD, FORCE) - total_cores_in_use}}
                        ]
                    }
                },
                {
                    '$addFields': {
                        'priority': {
                            '$cond': {
                                'if': {'$in': ['$mixing_proportion', [0.25, 0.5, 0.75]]},
                                'then': 1,  # Assign priority 1 to preferred values
                                'else': 2   # Assign priority 2 to other values
                            }
                        }
                    }
                },
                {
                    '$sort': {
                        'priority': 1,  # Prioritize by the custom field
                        'mixing_proportion': 1  # Then sort by mixing_proportion within the same priority
                    }
                },
                {
                    '$limit': 1  # Optional: Get only one document
                }
            ]))
                    # Check if a task was found
            if task:
                vasp_input()
                task = task[0]
                NOTASK = False
                task_id, folder_string, cores = task['_id'], task["folder_str"], task["cores"]
                target_directory = os.path.join(directory_vasp_input, folder_string)
                os.chdir(target_directory)
                # Start subprocess
                # Step 1: Make the script executable
                chmod_cmd = f"chmod +x ./script.sh"
                subprocess.run(chmod_cmd, shell=True, check=True)  # `check=True` ensures an exception is raised on failure
                time.sleep(0.1)
                # Step 2: Execute the script
                cmd = f"./script.sh"  # Since the script is now executable, you can directly call it
                process = subprocess.Popen(cmd, shell=True)
                processes[process.pid] = process # store running process in dictionary entry with pid
                cores_list[process.pid] = cores
                # Update MongoDB
                collection_calculation.update_one({'_id': task['_id']}, {'$set': {'RUNNING_pid': process.pid,'RUNNING_path': target_directory}, '$unset': {'TO_RUN': ""}})
                os.chdir("../../") # change directory to save .out file in correct directory
            else:
                if NOTASK == False:
                    print("No task found.")
                    vasp_input()
                    #create_mcsqs()
                    calculate_energies('VASP', 0.5)
                    NOTASK = True
                time.sleep(60)

        # Check status of running processes
        for pid, process in list(processes.items()):
            if process.poll() is not None:  # Process has finished
                task = collection_calculation.find_one({'RUNNING_pid': pid})
                if process.poll() == 0:
                    # Move the folder
                    try:
                        shutil.move(task["RUNNING_path"], directory_vasp_results)
                    except FileNotFoundError as e:
                        print(f"Error: {e}")
                    except PermissionError as e:
                        print(f"Permission denied: {e}")
                    # process result
                    process_output()
                    # Update MongoDB
                    collection_calculation.update_one(
                        {'RUNNING_pid': pid}, 
                        {
                            '$unset': {'RUNNING_pid': "", 'RUNNING_path': ""}
                        })
                    # Remove from tracking
                    del processes[pid]
                    del cores_list[pid]
                elif process.poll() == 2: # no convergence after 100 scf cycles
                    input("2")
                    del processes[pid]
                    del cores_list[pid]
                elif process.poll() == 3: # more than 40 bfgs cycles --> start with last CELL_PARAMETERS
                    input("3")
                    del processes[pid]
                    del cores_list[pid]
                    
                else: # Run was terminated
                    # TODO if process.poll() between 100 and 300: raise k_points? smearing?
                    print(process.poll())
                    input()
                    del processes[pid]
                    del cores_list[pid]
    except Exception as e:
        error_time = time.strftime('%Y-%m-%d %H:%M:%S')
        error_message = f"An error occurred at {error_time}:\n{str(e)}\n{traceback.format_exc()}\n"

        # Log the error to error.log, appending to the existing content
        with open("error.log", "a") as error_log:
            error_log.write(error_message)

        print(f"An error occurred: {e}")
        # Ensure that the processes array is still saved, and the while loop continues

    # Add a small delay to avoid a tight loop consuming excessive resources
    time.sleep(5)