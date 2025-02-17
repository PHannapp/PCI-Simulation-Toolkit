# Helper function to safely read a file and return its content or None
def safe_read_file(filepath):
    try:
        with open(filepath, "r") as file:
            return file.read()
    except FileNotFoundError:
        return None

def safe_read_lines(filepath):
    try:
        with open(filepath, "r") as file:
            return file.readlines()
    except FileNotFoundError:
        return None