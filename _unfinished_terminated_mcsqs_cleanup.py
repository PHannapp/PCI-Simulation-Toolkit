from pymongo import MongoClient
from bson.objectid import ObjectId  # Required to work with _id

# Connect to MongoDB
client = MongoClient('mongodb://localhost:27017/')
db = client.DFTforCALPHAD
collection_system = db.System
collection_calculation = db.calculation
input("Really do this?")
tasks = collection_calculation.find(
    {
        '$or': [
            {'RUNNING_pid': {'$exists': True}},
            {'FILES.POSCAR': {'$exists': True, '$eq': None}}  # Exists but is null
        ]
    }
)
input("Very sure?")
for task in tasks:
    updates = {}  # Dictionary to hold the update commands for MongoDB

    if task.get('RUNNING_pid'):
        # If RUNNING_pid exists, unset RUNNING_pid and RUNNING_path fields
        updates['$unset'] = {
            'RUNNING_pid': "",
            'RUNNING_path': ""
        }
        updates['$set'] = {'TO_RUN': True}  # Set TO_RUN to True

    if task.get('FILES'):
        # If FILES.POSCAR exists (even if it is None), unset FILES and set TO_RUN
        updates['$unset'] = {'FILES': ""}  # Unset the entire FILES field
        updates['$set'] = {'TO_RUN': True}  # Set TO_RUN to True
    # Apply the update to the document
    collection_calculation.update_one(
        {'_id': task['_id']},
        updates
    )