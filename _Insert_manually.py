from pymongo import MongoClient
import os
import shutil
from VASP_input import vasp_input

# Connect to MongoDB
client = MongoClient('mongodb://localhost:27017/')
db = client.DFTforCALPHAD
collection_system = db.System
collection_calculation = db.calculation

dir = "/home/pet69516/Dokumente/DFTforCALPHAD_home/POSCARS/_manually"
target_dir = "/home/pet69516/Dokumente/DFTforCALPHAD_home/POSCARS/_manually/old"

for file_path in os.listdir(dir):
    file_name = file_path.split('.')[0]
    POSCAR_path = os.path.join(dir, file_path)
    if os.path.isdir(POSCAR_path):
        continue
    with open(POSCAR_path, 'r') as file:
        file_content = file.read()
    lines = file_content.splitlines()
    elements = lines[5].split()
    numbers = lines[6].split()
    number_of_atoms = sum([int(n) for n in numbers])

    # Example document with POSCAR data
    task_data = {
        "name": file_name,
        "type": "manual",
        "POSCAR": file_content,
        "elements": elements,
        "number_of_atoms": number_of_atoms,
        "cores": 8,
        "TO_RUN": False
    }
    # Inserting the document
    collection_calculation.insert_one(task_data)
    # Move the file
    shutil.move(os.path.join(dir, file_path), os.path.join(target_dir, file_path))

vasp_input()