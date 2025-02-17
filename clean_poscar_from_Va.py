from pymongo import MongoClient
from bson.objectid import ObjectId  # Required to work with _id
from write import clean_mcsqs_poscar
import re

# Connect to MongoDB
client = MongoClient('mongodb://localhost:27017/')
db = client.DFTforCALPHAD
collection_system = db.System
collection_calculation = db.calculation


def insert_space_if_pattern_found(s):
    # Regex pattern to match "{int}-{int}"
    pattern = r'(\d+)-(\d+)'

    # If the pattern is found, insert space before the dash
    updated_string = re.sub(pattern, r'\1 -\2', s)
    
    return updated_string

def replace_xxx(s):
    if 'xxx' in s:
        return s.replace('xxx', '1.0')

tasks = collection_calculation.find({'POSCAR': {'$exists':True}})


def clean_poscar_from_Va():
    CONT = True
    while CONT:
        CONT = False
        tasks = collection_calculation.find(
            {
                "$or": [
                    {'POSCAR': {'$exists': True}},         # Condition 1: 'POSCAR' exists
                    {'FILES.POSCAR': {'$exists': True}}    # Condition 2: 'FILES.POSCAR' exists
                ]
            }
        )
        for task in tasks:
            if task.get("POSCAR"):
                uncleaned = task["POSCAR"]
            elif task.get("FILES") and task["FILES"].get("POSCAR"):
                uncleaned = task["FILES"]["POSCAR"]
            else:
                continue
            FOUND, poscar, number_of_atoms = clean_mcsqs_poscar(uncleaned)
            if FOUND:
                CONT = True
            collection_calculation.update_one(
                {'_id':ObjectId(task["_id"])}, 
                {
                    '$set': {'POSCAR': poscar,
                             'number_of_atoms': number_of_atoms}
                })
        

    tasks = collection_calculation.find(
        {
            "$or": [
                {'POSCAR': {'$exists': True}},         # Condition 1: 'POSCAR' exists
                {'FILES.POSCAR': {'$exists': True}}    # Condition 2: 'FILES.POSCAR' exists
            ]
        }
    )
    for task in tasks:
        if task.get('POSCAR'):
            poscar = insert_space_if_pattern_found(task['POSCAR'])
            collection_calculation.update_one(
            {'_id':ObjectId(task["_id"])}, 
            {
                '$set': {'POSCAR': poscar}
            })
