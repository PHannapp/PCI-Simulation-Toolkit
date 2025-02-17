import json
import time
import subprocess
import time
from pymongo import MongoClient, DESCENDING, ASCENDING
import shutil
import os
import traceback
from write import write_mcsqs, write_rndstr, turn_TO_RUN_false_ith
from read import safe_read_file, safe_read_lines
import psutil
from create_database_entries_from_system import create_reference_elements, create_endmembers, create_mcsqs, create_reference_elements_all
from clean_poscar_from_Va import clean_poscar_from_Va
from VASP_input import vasp_input, process_output, process_output_mcsqs, vasp_input_advanced, vasp_input_tmp
from CHGNet_helpers import save_trajectory, CHGNet_calculation

# Connect to MongoDB
client = MongoClient('mongodb://localhost:27017/')
db = client.DFTforCALPHAD
collection_system = db.System
collection_calculation = db.calculation

# Add a small delay to avoid a tight loop consuming excessive resources
#create_reference_elements_all()
#create_reference_elements()
#create_endmembers()
#create_mcsqs()
#clean_poscar_from_Va()
vasp_input()
#vasp_input_advanced('eff', 'HSE06', '671666f67a943236245787e0', "LHFCALC = .TRUE.\nGGA = PE\nHFSCREEN = 0.2")
#vasp_input_advanced('eff', 'GGA=91', '671666f67a943236245787e0')
#write_rndstr()
#clean_poscar_from_Va()
#process_output()
#vasp_input_tmp()
#process_output_mcsqs()
#turn_TO_RUN_false_ith(2)
#CHGNet_calculation(db.calculation_24_10_23, source='VASP.eff.', file="CONTCAR_SAVE2", type="manual") 
#CHGNet_calculation(db.calculation_24_10_23, source='VASP.eff.', file="CONTCAR_EFF", type="manual") 
#CHGNet_calculation(db.calculation, source='VASP.eff.', file="CONTCAR_EFF", type="manual") 
#CHGNet_calculation(db.calculation, source='VASP.eff.', file="CONTCAR_EFF") 
#CHGNet_calculation(db.calculation, source='VASP.prec.', file="CONTCAR_PREC") 
#CHGNet_calculation(db.calculation) 