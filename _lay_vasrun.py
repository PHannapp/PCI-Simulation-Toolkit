from pymongo import MongoClient
from bson.objectid import ObjectId  # Required to work with _id
import matplotlib.pyplot as plt
from pymatgen.core import Structure
import py3Dmol
from calculate import fit_interaction_parameters
import os

# Connect to MongoDB
client = MongoClient('mongodb://localhost:27017/')
db = client.DFTforCALPHAD
collection_system = db.System
collection_calculation = db.calculation


task_dir = r"/home/pet69516/Dokumente/DFTforCALPHAD/Phonopy"
def lay_vasprun(keyword1, keyword2, name = None):
    tasks = collection_calculation.find(
        {
            'type': "manual",
            keyword1: {'$exists': True}
        }
    )
    for task in tasks:
        if task.get('name'):
            if name not in task['name']:
                continue
            with open(os.path.join(task_dir, f"{task['name']}"), 'w') as POTCARfile:
                POTCARfile.write(task[keyword1][keyword2]['VASPRUN']['XML'])

lay_vasprun('VASP', 'general', name='Ni_')