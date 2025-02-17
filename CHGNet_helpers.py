import numpy as np
from pymatgen.core import Structure
from chgnet.model import CHGNet
from urllib.request import urlopen
from chgnet.model import StructOptimizer
from pymongo import MongoClient
from bson.objectid import ObjectId  # Required to work with _id
import pickle
from analyze import get_nested_value


def save_trajectory(trajectory):
    # Serialize the trajectory object into a binary format
    trajectory_data = {
        "energy": trajectory.energies,
        "forces": trajectory.forces,
        "stresses": trajectory.stresses,
        "magmoms": trajectory.magmoms,
        "atom_positions": trajectory.atom_positions,
        "cell": trajectory.cells,
        "atomic_number": trajectory.atoms.get_atomic_numbers(),
    }
    serialized_trajectory = pickle.dumps(trajectory_data)
    return serialized_trajectory


def CHGNet_calculation(collection_calculation, source="", file="POSCAR", type={'$exists': True}):
    tasks = collection_calculation.find(
        {
        f'{source}{file}': {'$exists': True, "$ne": None},
        'type': type,
        f'{source}CHGNet_results': {'$exists': False},
        }
    )
    list_tasks = list(tasks)
    len_tasks = len(list_tasks)
    print(f"Number of tasks: ", len_tasks)
    for ith_task ,task in enumerate(list_tasks):
        print(f"Task {ith_task} / {len_tasks}")
        print(task["elements"])
        print(task["_id"])
        if source != "":
            source_parts = source.split(".")[:-1]
            nested_source_value = get_nested_value(task, source_parts)
            structure = Structure.from_str(nested_source_value[file], fmt="poscar")
        else:
            structure = Structure.from_str(task[file], fmt="poscar")
        np.set_printoptions(precision=4, suppress=True)
        chgnet = CHGNet.load()
        relaxer = StructOptimizer()
        # Relax the perturbed structure
        result = relaxer.relax(structure, verbose=False)
        try:
            prediction = chgnet.predict_structure(structure)
            CHGNet_results = {
                # TODO automated determination of optimal DFT parameters (with ML?)
                "trajectory": save_trajectory(result["trajectory"]),
                "final_structure": result["final_structure"].to(fmt="poscar"),
                "type": "reference",
                "energy": prediction["e"].item(),
                "forces": prediction["f"].tolist(),
                "stress": prediction["s"].tolist(),
                "magnetic_moments": prediction["m"].tolist()
            }
        except ValueError as e:
            CHGNet_results = {
                # TODO automated determination of optimal DFT parameters (with ML?)
                "ERROR": str(e)
            }
            print("ERROR: ", e)
        collection_calculation.update_one(
        {'_id': ObjectId(task["_id"])},
        {
            '$set': {f'{source}CHGNet_results': CHGNet_results}
        })