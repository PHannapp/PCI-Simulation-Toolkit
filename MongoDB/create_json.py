from connect import collection_calculation, collection_system, db
import copy

collection_system_new = db.collection_system_pub
collection_calculation_new = db.collection_calculation_pub

systems_elements = {
    "BCC": [[["NB"], ["H"]], [["TI"], ["H"]], [["VV"], ["H"]]],
    "FCC": [[["NB"], ["H"]], [["TI"], ["H"]], [["VV"], ["H"]]],
    "C14": [[["TI"], ["MN"], ["H"]]],
    "AB5": [[["LA"], ["NI"], ["H"]]],
    "ABalpha": [[["TI"], ["FE"], ["H"]]]
}

find_set_tmp = {"extras": {"$exists": False}, "timestamp": {"$exists": False}, "type": "endmember"}

fields_to_delete = ["CHGNet", "CHGNet_contcar", "folder_str", "VASP.POTCAR", "VASP.SCRIPT"]

def copy_entries():
    for key, element_sets in systems_elements.items():
        print(key)
        # Find and copy system document
        task_system = collection_system.find_one({"system": key})
        if task_system:
            # Remove _id to avoid duplication issues
            task_system_copy = copy.deepcopy(task_system)
            task_system_copy.pop("_id", None)
            collection_system_new.insert_one(task_system_copy)

        # Process each element combination
        for elements in element_sets:
            print(elements)
            find_set = copy.deepcopy(find_set_tmp)
            find_set["system_name"] = key
            find_set["elements"] = elements

            tasks = collection_calculation.find(find_set)
            for task in tasks:
                task_copy = copy.deepcopy(task)
                task_copy.pop("_id", None)
                collection_calculation_new.insert_one(task_copy)
    
    tasks = collection_calculation.find({"type": "reference"})
    for task in tasks:
        task_copy = copy.deepcopy(task)
        task_copy.pop("_id", None)
        collection_calculation_new.insert_one(task_copy)

def clean():
    tasks = collection_calculation_new.find({})
    for task in tasks:
        unset_fields = {}
        for field in fields_to_delete:
            # Check if the field exists in the document
            if field in task:
                unset_fields[field] = ""  # Value doesn't matter, presence of key triggers removal
        if unset_fields:
            collection_calculation_new.update_one(
                {"_id": task["_id"]},
                {"$unset": unset_fields}
            )

#copy_entries()
clean()

