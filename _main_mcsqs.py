import json
import time
import subprocess
import time
from pymongo import MongoClient, DESCENDING, ASCENDING
import shutil
import os
import traceback
from write import write_mcsqs, write_rndstr
from read import safe_read_file, safe_read_lines
import psutil
from create_database_entries_from_system import create_reference_elements, create_endmembers, create_mcsqs
from clean_poscar_from_Va import clean_poscar_from_Va
from VASP_input import vasp_input, process_output

# Save the current working directory
original_directory = os.getcwd()

def kill_process_and_children(pid):
    try:
        parent = psutil.Process(pid)
        children = parent.children(recursive=True)  # Get child processes
        for child in children:
            child.kill()  # Kill child processes
        parent.kill()  # Kill the parent process
        parent.wait()  # Wait for the parent to be fully terminated
    except psutil.NoSuchProcess:
        print(f"Process {pid} does not exist")

# 1. 

#convergence(testfile, calculation, 'Ce')
def read_config():
    default_values = {'THRESHOLD': 0, 'FORCE': 1000}
    try:
        with open(os.path.join(original_directory, 'config_mcsqs.json'), 'r') as file:
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

def process_result(collection, pid):
    task = collection.find_one({'RUNNING_pid': pid})
    # TODO copy outputfile in MongoDB and delete everything
    with open(f"{task['RUNNING_path']}/{task['name']}.out", "r") as file:
        file_contents = file.read()
        collection.update_one({'_id': task['_id']}, {'$set': {'OUTPUT_FILE': file_contents}})
    #shutil.rmtree(task['RUNNING_path'])
    task = collection.find_one({'RUNNING_pid': pid})
    #if task['cal']['type'] == 'formation_enthalpy':
    lines_cv, lines_ap, k_point_value, energy, convergence, cpu_time, wall_time, cell_volume, density = extract_vc(task['OUTPUT_FILE'])
    collection.update_one(
        {'_id': task['_id']},
        {
            '$set': {
                'CELL_PARAMETERS_ALAT': lines_cv,
                'ATOMIC_POSITIONS_crystal': lines_ap,
                'K_POINT_VALUE': k_point_value,
                'ENERGY': energy,
                'CONVERGENCE': convergence,
                'CPU_TIME': cpu_time,
                'WALL_TIME': wall_time,
                'CELL_VOLUME': cell_volume,
                'DENSITY': density
            }
        }
    )
    
# Connect to MongoDB
client = MongoClient('mongodb://localhost:27017/')
db = client.DFTforCALPHAD
collection_system = db.System
collection_calculation = db.calculation

processes = {} # dictionary to store running processes
cores_list = {} # total cores in usage
NOTASK = False
t_save = time.time()
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
            collection_calculation.update_one({'RUNNING_pid': pid_to_terminate}, {'$set': {'TO_RUN': False}, '$unset': {'RUNNING_pid': "", 'RUNNING_path': "", 'INPUT_FILE': ""}})
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
                            {'TO_RUN': False},
                            {'type': 'mcsqs'},
                            {'RNDSTR': {'$exists': True}}
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
                task = task[0]
                NOTASK = False

                # Start subprocess
                #cmd = create_command(task_id)  # Function to create the command based on the task
                full_file_path = os.path.join('./RUNS/', task["type"] + "_" + str(task['_id']))
                if not os.path.exists(full_file_path):
                    os.makedirs(full_file_path)
                os.chdir(full_file_path) # change directory to save .out file in correct directory
                # write RNDSTR.in file
                with open(f"./rndstr.in", "w") as f:
                    f.write(task ["RNDSTR"])
                with open(f"./sqscell.out", "w") as f:
                    f.write(task ["SQSCELL"])
                with open(f"./{str(task['_id'])}.py", "w") as f:
                    f.write(write_mcsqs(task))
                env = os.environ.copy()  # Copy current environment variables
                process = subprocess.Popen(f"python3 {str(task['_id'])}.py", shell=True, env=env) # process containts information like pid
                #time.sleep(2)
                #child_pid = subprocess.check_output(f"pgrep -P {process.pid}", shell=True).decode().strip()
                #if child_pid.isdigit():
                #    processes[int(child_pid)] = process
                processes[process.pid] = process # store running process in dictionary entry with pid
                cores_list[process.pid] = 1# task["calculation_information"]["cores"]
                # Update MongoDB
                collection_calculation.update_one({'_id': task['_id']}, {'$set': {'RUNNING_pid': process.pid,'RUNNING_path': full_file_path}, '$unset': {'TO_RUN': ""}})
                os.chdir("../../") # change directory to save .out file in correct directory
            else:
                if NOTASK == False:
                    print("No task found.")
                    NOTASK = True
                clean_poscar_from_Va()
                time.sleep(60)

        # Check status of running processes
        for pid, process in list(processes.items()):
            if process.poll() is not None:  # Process has finished
                task = collection_calculation.find_one({'RUNNING_pid': pid})
                if process.poll() == 0: # Run ended gracefully
                    # process result
                    #process_result(collection_calculation, pid)
                    # Read and print the contents of stdout.txt
                    # Read files safely
                    # Handle the specific POSCAR file case
                    CONTINUE = False
                    for i in range(1):
                        lines = safe_read_lines(os.path.join(task['RUNNING_path'], "bestsqs.out-POSCAR"))
                        if lines is not None and len(lines) >= 2:
                            lines[1] = "1.0\n"  # Replace the second line (index 1) with "1.0"
                            POSCAR = ''.join(lines)
                            CONTINUE = True
                            break
                        else:
                            print(i, "th time no bestsqs.out-POSCAR found!")
                            time.sleep(10)
                    if not CONTINUE:
                        POSCAR = None

                    stdout_content = safe_read_file(os.path.join(task['RUNNING_path'], "stdout.txt"))
                    stderr_content = safe_read_file(os.path.join(task['RUNNING_path'], "stderr.txt"))
                    stdout_2_content = safe_read_file(os.path.join(task['RUNNING_path'], "stdout_main.txt"))
                    stderr_2_content = safe_read_file(os.path.join(task['RUNNING_path'], "stderr_main.txt"))
                    bestcorr = safe_read_file(os.path.join(task['RUNNING_path'], "bestcorr.out"))
                    bestsqs = safe_read_file(os.path.join(task['RUNNING_path'], "bestsqs.out"))
                    clusters = safe_read_file(os.path.join(task['RUNNING_path'], "clusters.out"))
                    mcsqs_log = safe_read_file(os.path.join(task['RUNNING_path'], "mcsqs.log"))
                    rndstrgrp = safe_read_file(os.path.join(task['RUNNING_path'], "rndstrgrp.out"))
                    sqscell = safe_read_file(os.path.join(task['RUNNING_path'], "sqscell.out"))

                    # Create the dictionary with the safely read files
                    files = {
                        'STDOUT': [stdout_content, stdout_2_content],
                        'STDERR': [stderr_content, stderr_2_content],
                        'bestcorr': bestcorr,
                        'bestsqs': bestsqs,
                        'clusters': clusters,
                        'mcsqs_log': mcsqs_log,
                        'rndstrgrp': rndstrgrp,
                        'sqscell': sqscell,
                        'POSCAR': POSCAR
                    }
                    # Update MongoDB
                    print("Saved task with task_id:", task['_id'])
                    collection_calculation.update_one(
                        {'RUNNING_pid': pid}, 
                        {
                            '$unset': {'RUNNING_pid': "", 'RUNNING_path': ""},
                            '$set': {'FILES': files,
                                     'TO_RUN': False}
                        })
                    
                    vasp_input()
                    # Remove from tracking    
                    # # Kill any remaining processes (parent and children)
                    kill_process_and_children(pid)
                    
                    del processes[pid]
                    del cores_list[pid]
                else: # Run was terminated
                    collection_calculation.update_one(
                        {'RUNNING_pid': pid},
                        {
                            '$unset': {'RUNNING_pid': "", 'RUNNING_path': ""},
                            '$set': {'TERMINATED': process.poll()}
                        })
                    del processes[pid]
                    del cores_list[pid]
        #if time.time() - t_save > 60*60*2:
        #    t_save = time.time()
        #    create_reference_elements()
        #    create_endmembers()
        #    create_mcsqs()
        #    clean_poscar_from_Va()
        #    vasp_input()
        #    write_rndstr()
        #    process_output()
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