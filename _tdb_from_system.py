from pymongo import MongoClient
from bson.objectid import ObjectId  # Required to work with _id
import os
from analyze import get_nested_value

# Connect to MongoDB
client = MongoClient('mongodb://localhost:27017/')
db = client.DFTforCALPHAD
collection_system = db.System
collection_calculation = db.calculation

calculation_collection_name = collection_calculation.name

system_id = '678a29768c3e8eaa28454296'
sources = ["VASP"]#METAGGA=R2SCAN.eff", "VASP.prec", 'VASP.METAGGA=R2SCAN.prec']#, 'VASP.METAGGA=R2SCAN.prec', "VASP.GGA=91.eff", "VASP.GGA=91.prec"]
#"VASP.eff", "VASP.prec", 
for source in sources:
    task_system = collection_system.find_one(
        {
            '_id': ObjectId(system_id)
        }
    )
    database_name = ""
    tmp_str = ""
    unique_elements = []
    ## create the header string and a unique list of all different elements
    for s, sublattice in enumerate(task_system['elements']):
        for i, element in enumerate(sublattice):
            if i != 0:
                tmp_str += '-'
            if element == 'H,Va':
                tmp_str += 'H-Va'
                unique_elements.append('H2')
            else:
                tmp_str += element
                unique_elements.append(element)
        tmp_str += "_"
        tmp_str += str(task_system['multiplicity'][s])
        database_name += tmp_str
        tmp_str = "_"
    database_name += "_" + source.split('_')[0]
    unique_elements = list(set([element.upper() for element in unique_elements]))
    # 1. create file and name, import numpy
    final_string = "import numpy as np\n# reference -------------------------------------------------------------------------------\n"
    # 2. create ref
    for element in unique_elements:
        with open(os.path.join("CALPHAD_FILES", 'CONSTANTS', 'GHSER.py'), 'r') as file:
            for line in file:
                if f'GHSER{element}' in line:
                    # Append the line
                    final_string += line
                    # Try to get the next line, if it exists
                    next_line = next(file, None)
                    if next_line:
                        final_string += next_line

    final_string += "\n# Endmembers -------------------------------------------------------------------------------\n"
    # 3. create endmembers
    tasks = collection_calculation.find(
        {
            'system_id': ObjectId(system_id),
            'type': 'endmember',
            f'{source}.formation_enthalpy_for_calphad': {'$exists': True}
        }
    )

    # Aggregation pipeline
    pipeline = [
        {
            # Match documents based on the specified conditions
            "$match": {
                'system_id': ObjectId(system_id),
                "type": "endmember",
                f"{source}.formation_enthalpy_for_calphad": {"$exists": True}
            }
        },
        {
            "$group": {
                "_id": {
                    "elements": "$elements",
                    "mixing_proportion": "$mixing_proportion"
                },
                "min_eff_energy": {"$min": f"${source}.formation_enthalpy_for_calphad"},
                "doc_id": {"$first": "$_id"}  # Capture the ID of one document with the lowest energy
            }
        },
        {
            # Debug: Add a stage to output intermediate results
            "$project": {
                "elements": "$_id.elements",
                "mixing_proportion": "$_id.mixing_proportion",
                "min_eff_energy": 1,
                "doc_id": 1
            }
        },
        {
            "$lookup": {
                "from": calculation_collection_name,  # Replace with the actual collection name
                "localField": "doc_id",
                "foreignField": "_id",
                "as": "lowest_energy_doc"
            }
        },
        {
            "$unwind": {
                "path": "$lowest_energy_doc",
                "preserveNullAndEmptyArrays": True  # Set to True to identify mismatches
            }
        },
        {
            "$replaceRoot": {
                "newRoot": "$lowest_energy_doc"
            }
        }
    ]

    # Execute the pipeline
    results = list(collection_calculation.aggregate(pipeline))
    #def GCE_NI_V_V(T):
    #    return -197460 + (24.3084) * T + GHSERCE(T) + 5*GHSERNI(T)
    for task in results:
        line1 = 'def G'
        line2 = '   return '
        # for the calculation of the entropy
        amount_of_Hs = 0
        for sublattice, element in enumerate(task['elements']):
            if sublattice != 0:
                line1 += '_'
            if element[0] == 'H':
                line1 += 'H'
                amount_of_Hs += task_system['multiplicity'][sublattice] / 2
            elif element[0] == 'Va':
                line1 += 'V'
            else:
                line1 += element[0].upper()
                line2 += f'{task_system["multiplicity"][sublattice]}*GHSER{element[0].upper()}(T) + '
        line2 += f'{amount_of_Hs}*GHSERH2(T) + '
        hydrogen_entropy = amount_of_Hs * 130

        source_parts = source.split(".")
        nested_source_value = get_nested_value(task, source_parts)
        formation_enthalpy = round(nested_source_value['formation_enthalpy_for_calphad'] * 1e3)
        line2 += f'({formation_enthalpy}) + ({hydrogen_entropy} * T)\n'
        line1 += '(T):\n'
        final_string += line1 + line2

    final_string += "\n# Interaction parameters -------------------------------------------------------------------\n"

    # 5. write interaction parameters in tdb

    tasks = collection_calculation.find(
        {
            'system_id': ObjectId(system_id),
            'type': 'mcsqs',
            f'{source}.interaction_parameters': {'$exists': True}
        }
    )
    #def LCE_NI_HV_H_0(T):
    #    return -20632
    task_elements = []
    for task in tasks:
        source_parts = source.split(".")
        nested_source_value = get_nested_value(task, source_parts)
        if not nested_source_value.get("number_of_atoms_per_formula_unit"):
            continue
        if task_elements != task['elements']:
            task_elements = task['elements']
        else:
            continue
        line1 = 'def L'
        # for the calculation of the entropy
        for sublattice, element in enumerate(task['elements']):
            if sublattice != 0:
                line1 += '_'
            if element == ['H', 'Va']:
                line1 += 'HV'
            elif len(element) == 2:
                line1 += element[0].upper()
                line1 += element[1].upper()
            elif element[0] == 'Va':
                line1 += 'V'
            else:
                line1 += element[0].upper()
        for order, interaction_parameter in enumerate(nested_source_value["interaction_parameters"]):
            line1_tmp = line1 + f'_{order}(T):\n'
            line2 = f'   return {interaction_parameter }\n'#/ nested_source_value["number_of_atoms_per_formula_unit"]}\n'
            final_string += line1_tmp + line2
    #def LCE_NI_HV_H_0(T):
    #    return -20632


    #### NOT YET
    # 6. calculate Entropy

    collection_system.update_one(
        {
            '_id': ObjectId(system_id)
        },
        {
            '$set': {f'{source}_tdb': final_string} 
        }
    )
    source_parts = source.split('.')
    source_db = '_'.join(source_parts)
    with open(os.path.join("CALPHAD_CALC", 'database.py'), 'w') as file:
        file.write(final_string)
    with open(os.path.join("CALPHAD_CALC", f'database_{source_db}.py'), 'w') as file:
        file.write(final_string)
    print(final_string)