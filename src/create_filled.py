from src.read import read_poscar
from src.calculate_NN import seperate_vac_atoms
import os
import numpy as np
from itertools import combinations, product
import random
from pymongo import MongoClient

def create_NN():
    # Connect to MongoDB
    client = MongoClient('mongodb://localhost:27017/')
    db = client.DFTforCALPHAD
    collection_system = db.System
    collection_calculation = db.calculation

    tasks = collection_calculation.find()
    vacs = {}
    for task in tasks:
        if task["type"] == "manual":
            name_array = task['name'].split("_")
            if len(name_array) == 7:
                mat_class = ("_").join(name_array[1:3])
            else:
                mat_class = name_array[1]
            vac = ("_").join(name_array[-2:])
            if mat_class not in vacs:
                vacs[mat_class] = {}  # Create new list
            if vac in vacs[mat_class]:
                vacs[mat_class][vac].append(task['VASP']['energy'])  # Extend existing list
            else:
                vacs[mat_class][vac] = [task['VASP']['energy']]  # Create new list
    for mat_class, vac in vacs.items():
        for vac_pos, energies in vac.items():
            print(f"{mat_class}'s {vac_pos} interstitital has energies: {energies}\n")
            #contcar = task["VASP.CONTCAR"]
            #oszicar = task["VASP.OSZICAR"]



def update_NN():
    lattice_vectors, atoms = read_poscar(VACS=True, file_type = 'NN')
    for atom in atoms:
        print(atom.label)
        print(atom.position)    
    vac_atoms = seperate_vac_atoms(atoms)
    for vac_atom in vac_atoms:
        vac_atom.position = np.array([float(pos) for pos in vac_atom.position])
    # Generate the translation vectors using itertools.product
    translation_vectors = np.array(list(product(np.arange(0, 1), np.arange(0, 1), np.arange(0, 1))))
    lowestt_distance = [1e9, 0, 0]
    for ith_atom in range(len(vac_atoms)):
        for jth_atom in range(ith_atom + 1, len(vac_atoms)):
            lowest_distance = 1e9
            for translation_vector in translation_vectors:
                vac_atom1 = np.dot((vac_atoms[ith_atom].position + translation_vector), lattice_vectors)
                vac_atom2 = np.dot((vac_atoms[jth_atom].position + translation_vector), lattice_vectors)
                distance = np.linalg.norm(vac_atom1 - vac_atom2)
                if distance < lowest_distance:
                    lowest_distance = distance
            print(ith_atom, jth_atom, lowest_distance)
            if lowest_distance < lowestt_distance[0]:
                lowestt_distance[0] = lowest_distance
                lowestt_distance[1] = ith_atom
                lowestt_distance[2] = jth_atom
    print(lowestt_distance[1], lowestt_distance[2], lowestt_distance[0])
            
