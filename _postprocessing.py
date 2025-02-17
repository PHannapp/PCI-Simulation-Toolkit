from analyze import print_energies, print_energies_from_OSZICAR, plot_endmember_energies, plot_mixing_energies, show_strcuture, fit_mixing_energies_to_interaction_parameters
from calculate import calculate_energies

# Save the plot to a directory with a specified name
output_dir = "/home/pet69516/Dokumente/DFTforCALPHAD_home/PLOTS"  # Specify your directory here
sources = ["VASP"]
force_thr = 0.5

#print_energies_from_OSZICAR({'$exists':False}, 'VASP', 'eff')
for source in sources:
    #calculate_energies(sources, force_thr)
    #input()
    #print_energies("(A)1,(B)2,(H,Va)3-1/3,(H,Va)3-1/3,(H,Va)6-1/3", 'VASP', 'eff', energy='formation_enthalpy')
    #input()
    plot_endmember_energies(source, 'formation_enthalpy', force_thr, 'VV', save=output_dir)
    ##
    #input()
    #
    #input()
    fit_mixing_energies_to_interaction_parameters("(A),(B)2(B)3,(H,Va)1,(H,Va)6", source, bound=4e6, force=force_thr, parameter_space=[1,0,0], save=output_dir)
    ##show_strcuture()
    plot_mixing_energies("(A),(B)2(B)3,(H,Va)1,(H,Va)6", source, 'mixing_enthalpy', force_thr, exclude='VV', save=output_dir)
    pass

