from src.VASP_input import vasp_input, process_output

# Add a small delay to avoid a tight loop consuming excessive resources
vasp_input()
process_output()