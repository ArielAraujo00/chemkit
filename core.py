# Built-in
import re
import json
from pathlib import Path

# --------------- constants --------------- #
_EH_TO = {'kcal/mol': 627.50961,
          'kJ/mol': 2625.5002,
          'eV': 27.211399,
          'cm-1': 219474.63,
          'hartree': 1}

_TIME_TO = {'s': 1,
            'ms': 1000,
            'min': 1/60,
            'h': 1/3600,
            'd': 1/86400}

# --------------- regex patterns --------------- #
_INT_RE   = re.compile(r'[+-]?\d+')

_FLOAT_RE = re.compile(r'[+-]?\d+\.\d+')

_SCY_RE   = re.compile(r'[+-]?\d*\.\d+[Ee][+-]?\d+')

_NUMBER_RE = re.compile(f'{_SCY_RE.pattern}|{_FLOAT_RE.pattern}|{_INT_RE.pattern}')


# --------------- main helpers --------------- #
def parse_path(file_path):
    assert isinstance(file_path, (str, Path)), 'Must be a path-like object (str | Path)'
    if not isinstance(file_path, Path):
        return Path(file_path)
    else:
        return file_path

def read_lines(file):
    file = parse_path(file)
    with file.open('r') as f:
        return f.readlines()

def extract_numbers(line):
    out = []
    for s in _NUMBER_RE.findall(line):
        if 'e' in s or 'E' in s:
            out.append(float(s))
        elif '.' in s:
            out.append(float(s))
        else:
            out.append(int(s))
    return out

def safe_float(number):
    try:
        return float(number)
    except (ValueError, TypeError):
        return None

def safe_int(number):
    try:
        return int(number)
    except (ValueError, TypeError):
        return None


# --------------- data storage --------------- #
def dump_json(path, data):
    path = Path(path)
    with path.open("w") as f:
        json.dump(data, f, indent=2)

def load_json(path):
    path = Path(path)
    with path.open("r") as f:
        return json.load(f)

def load_descriptor_dataframe(path):
    json_files = list(Path(path).glob('*.json'))
    def _rows():  # Build from generator -> no need to load all files into memory
        for file in tqdm(json_files, total=len(json_files)):
            data = load_json(file)
            data['__index__'] = file.stem
            yield data
    return pd.DataFrame(_rows()).set_index('__index__')

# --------------- refactor latter --------------- #
"""
import os, shutil, glob
from decimal import Decimal

################################################################################ Text editing
# function to convert to superscript
def super_script(string):
    normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
    super_s = "ᴬᴮᶜᴰᴱᶠᴳᴴᴵᴶᴷᴸᴹᴺᴼᴾQᴿˢᵀᵁⱽᵂˣʸᶻᵃᵇᶜᵈᵉᶠᵍʰᶦʲᵏˡᵐⁿᵒᵖ۹ʳˢᵗᵘᵛʷˣʸᶻ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾"
    res = string.maketrans(''.join(normal), ''.join(super_s))
    return string.translate(res)

# function to convert to subscript
def sub_script(string):
    normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
    sub_s = "ₐ₈CDₑբGₕᵢⱼₖₗₘₙₒₚQᵣₛₜᵤᵥwₓᵧZₐ♭꜀ᑯₑբ₉ₕᵢⱼₖₗₘₙₒₚ૧ᵣₛₜᵤᵥwₓᵧ₂₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎"
    res = string.maketrans(''.join(normal), ''.join(sub_s))
    return string.translate(res)

# function to format the xyz line upon file conversion
def format_xyz(atom, coord, decimal=5):
    for i in range(len(coord)):
        if coord[i] < 0:
            coord[i] = f'{coord[i]:.{decimal}f}'
        else:
            coord[i] = f' {coord[i]:.{decimal}f}'
    if len(atom) == 1:
        atom = f' {atom} '
    else:
        atom = f' {atom}'
    return f'{atom}    {coord[0]}    {coord[1]}    {coord[2]}\n'

# makes a list of formated numbers
def format_range(start, end, step=1):
    out = list()
    # Format input
    start = Decimal(str(start))
    end = Decimal(str(end))
    step = Decimal(str(step))
    # Loop range
    current = start
    while current <= end:
        if current % 1 == 0:
            out.append(int(current))
        else:
            out.append(float(current))
        current += step
    return out

################################################################################ File & System

# function to retrieve files and/or directories paths from a given root
def get_main(path, select):
    assert select in ['dirs', 'files', 'all'], "select must be: \'dirs\', \'files\' or \'all\'"
    for root, dirs, files in os.walk(path):
        if root == path:
            if select == 'dirs':
                return [os.path.join(root, d) for d in dirs]
            if select == 'files':
                return [os.path.join(root, f) for f in files]
            if select == 'all':
                return [os.path.join(root, i) for i in os.listdir(root)]

# function to retrieve files and/or directories that match a given pattern string
def glob_files(path, pattern, select='all'):
    if select == 'all': # loop through all inner dirs and files of a root
        output = list()
        for root, dirs, files in os.walk(path):
            path_pattern = os.path.join(root, pattern)
            output.extend(glob.glob(path_pattern))
        return output
    elif select == 'root': # only the root path
        path_pattern = os.path.join(root, pattern)
        return glob.glob(path_pattern)

# makes sure the path exists and creates/clears it when don't
def secure_path(path, clear=False):
    if not os.path.exists(path):
        os.makedirs(path)
    elif os.path.exists(path) and clear:
        shutil.rmtree(path)
        os.makedirs(path)
    else:
        while True:
            ans = input(f'<{path}> already exists, do you want to overwirte? [y/n]')
            if ans == 'y':
                shutil.rmtree(path)
                os.makedirs(path)
                break
            elif ans == 'n':
                break
            else:
                print('Wrong input')

# extracts only file name from path string (removes extension)
def get_name(path):
    return os.path.splitext(os.path.split(path)[1])[0]

################################################################################ OTHERS
def energy_from_xyz(path):
    with open(path) as file:
        comment = file.readlines()[1]
        energy = re.findall(float_regex, comment)[0]
        return float(energy)
"""
