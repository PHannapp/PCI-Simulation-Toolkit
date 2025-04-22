from src.calculate import calculate_energies

sources = ["VASP"]
force_thr = 0.5

#print_energies_from_OSZICAR({'$exists':False}, 'VASP', 'eff')
for source in sources:
    calculate_energies(sources, force_thr)
    pass

