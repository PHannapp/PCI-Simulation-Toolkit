#import json
#import re
#
#import json
#
#with open("DFTforCALPHAD.collection_calculation_pub_3.json", "r") as f:
#    data = json.load(f)
#
#def print_keys_recursive(d, prefix=""):
#    if isinstance(d, dict):
#        for k, v in d.items():
#            full_key = f"{prefix}.{k}" if prefix else k
#            print(full_key)
#            print_keys_recursive(v, full_key)
#    elif isinstance(d, list):
#        for i, item in enumerate(d):
#            print_keys_recursive(item, f"{prefix}[{i}]")
#
## Nur die ersten paar Einträge anschauen
#for i, entry in enumerate(data[:5]):
#    if "VASP" in entry:
#        print(f"\n--- Keys in entry {i} (VASP) ---")
#        print_keys_recursive(entry["VASP"])
#
#input()

import json
import re

# Lade die Datei
with open("DFTforCALPHAD.collection_calculation_pub_3.json", "r") as f:
    data = json.load(f)

pattern = re.compile(r"^SLURM-\d+(\..+)?$")  # z. B. SLURM-12345678 oder SLURM-12345678.OUT

total_removed = 0
cleaned_data = []

for entry in data:
    if isinstance(entry, dict) and "VASP" in entry and isinstance(entry["VASP"], dict):
        vasp_block = entry["VASP"]
        keys_before = set(vasp_block.keys())

        # Entferne alle Keys, die dem Muster entsprechen
        entry["VASP"] = {k: v for k, v in vasp_block.items() if not pattern.match(k)}

        keys_after = set(entry["VASP"].keys())
        removed = keys_before - keys_after
        if removed:
            print(f"Removed: {removed}")
            total_removed += len(removed)

    cleaned_data.append(entry)

print(f"\nTotal SLURM* keys removed: {total_removed}")

# Gespeicherte Datei
with open("DFTforCALPHAD.collection_calculation_pub_3_cleaned.json", "w") as f:
    json.dump(cleaned_data, f, indent=2)

