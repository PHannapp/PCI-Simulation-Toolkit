from pymongo import MongoClient

# Step 1: Connect to MongoDB
client = MongoClient("mongodb://localhost:27017/")  # Replace with your MongoDB connection string if different
db = client["DFTforCALPHAD"]  # Replace with your database name

# Step 2: Choose source and destination collections
source_collection = db["calculation"]  # Replace with your source collection name
destination_collection = db["calculation_241002"]  # Replace with your destination collection name

# Helper function to replace dots in keys
def replace_dots_in_keys(document):
    new_document = {}
    for key, value in document.items():
        new_key = key.replace('.', '_')  # Replace dots with underscores (or any other valid character)
        if isinstance(value, dict):  # If the value is a nested document, apply recursively
            value = replace_dots_in_keys(value)
        new_document[new_key] = value
    return new_document

# Step 3: Copy all documents from source collection to destination collection
def copy_collection(source_collection, destination_collection):
    documents = source_collection.find()  # Fetch all documents
    # Check if there are any documents to copy
    if source_collection.count_documents({}) > 0:
        processed_documents = [replace_dots_in_keys(doc) for doc in documents]  # Process documents to remove dots
        destination_collection.insert_many(processed_documents)  # Insert them into the destination collection

# Step 4: Execute the copy for "calculation" collection
copy_collection(source_collection, destination_collection)

# Step 5: Repeat for "System" collection
source_collection = db["System"]  # Replace with your source collection name
destination_collection = db["System_241002"]  # Replace with your destination collection name
copy_collection(source_collection, destination_collection)

print("Documents copied from 'calculation' and 'System' to their respective new collections.")
