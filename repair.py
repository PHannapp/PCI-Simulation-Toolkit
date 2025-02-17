from pymongo import MongoClient

# Step 1: Connect to MongoDB
client = MongoClient("mongodb://localhost:27017/")  # Replace with your MongoDB connection string if different
db = client["DFTforCALPHAD"]  # Replace with your database name

destination_collection = db["calculation"]  # Replace with your destination collection name

# Step 3: Copy all documents from source collection to destination collection
def copy_collection(destination_collection):
    documents = destination_collection.find({"type": "endmember"})  # Fetch all documents
    for document in documents:
        if document["VASP"]["eff"]["POSCAR"]:
            destination_collection.update_one(
                {"_id": document["_id"]},
                {"$set": {f"POSCAR": document["VASP"]["eff"]["POSCAR"]}}
            )


# Step 4: Execute the copy
copy_collection(destination_collection)
