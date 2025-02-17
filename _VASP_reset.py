from pymongo import MongoClient
from bson.objectid import ObjectId  # Required to work with _id
import os
import numpy as np
from write import MAGNETIZATION_DICT

# Connect to MongoDB
client = MongoClient('mongodb://localhost:27017/')
db = client.DFTforCALPHAD
collection_system = db.System
collection_calculation = db.calculation


directory = "/home/pet69516/Dokumente/DFTforCALPHAD/VASPInput"
POTCAR_dir = "/home/pet69516/Dokumente/DFTforCALPHAD/POTCARS"




tasks = collection_calculation.find(
    {
        "$and": [
            {'VASP': {'$exists': True}},         # 
            {'type': 'mcsqs'}   #
        ]
    }
)
for task in tasks:
    collection_calculation.update_one(
        {"_id": task["_id"]},
        {
            "$set": {
                'TO_RUN':True}
            ,
            "$unset": {
                "VASP": ""
            }
        }
    )
