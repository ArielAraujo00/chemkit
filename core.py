# Built-in
import re
import json
from pathlib import Path

import numpy as np

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

_ALPHABET_LC = 'abcdefghijklmnopqrstuvwxyz'

_ALPHABET_UC = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

_SUPER_MAP = str.maketrans({
    'A':'ᴬ', 'B':'ᴮ', 'C':'ᶜ', 'D':'ᴰ', 'E':'ᴱ', 'F':'ᶠ', 'G':'ᴳ',
    'H':'ᴴ', 'I':'ᴵ', 'J':'ᴶ', 'K':'ᴷ', 'L':'ᴸ', 'M':'ᴹ', 'N':'ᴺ',
    'O':'ᴼ', 'P':'ᴾ', 'Q':'Q', 'R':'ᴿ', 'S':'ˢ', 'T':'ᵀ', 'U':'ᵁ',
    'V':'ⱽ', 'W':'ᵂ', 'X':'ˣ', 'Y':'ʸ', 'Z':'ᶻ',

    'a':'ᵃ', 'b':'ᵇ', 'c':'ᶜ', 'd':'ᵈ', 'e':'ᵉ', 'f':'ᶠ', 'g':'ᵍ',
    'h':'ʰ', 'i':'ᶦ', 'j':'ʲ', 'k':'ᵏ', 'l':'ˡ', 'm':'ᵐ', 'n':'ⁿ',
    'o':'ᵒ', 'p':'ᵖ', 'q':'۹', 'r':'ʳ', 's':'ˢ', 't':'ᵗ', 'u':'ᵘ',
    'v':'ᵛ', 'w':'ʷ', 'x':'ˣ', 'y':'ʸ', 'z':'ᶻ',

    '0':'⁰', '1':'¹', '2':'²', '3':'³', '4':'⁴',
    '5':'⁵', '6':'⁶', '7':'⁷', '8':'⁸', '9':'⁹',

    '+':'⁺', '-':'⁻', '=':'⁼', '(':'⁽', ')':'⁾'
})


_SUB_MAP = str.maketrans({
    'A':'ₐ', 'B':'₈', 'C':'C', 'D':'D', 'E':'ₑ', 'F':'բ', 'G':'G',
    'H':'ₕ', 'I':'ᵢ', 'J':'ⱼ', 'K':'ₖ', 'L':'ₗ', 'M':'ₘ', 'N':'ₙ',
    'O':'ₒ', 'P':'ₚ', 'Q':'Q', 'R':'ᵣ', 'S':'ₛ', 'T':'ₜ', 'U':'ᵤ',
    'V':'ᵥ', 'W':'w', 'X':'ₓ', 'Y':'ᵧ', 'Z':'Z',

    'a':'ₐ', 'b':'♭', 'c':'꜀', 'd':'ᑯ', 'e':'ₑ', 'f':'բ', 'g':'₉',
    'h':'ₕ', 'i':'ᵢ', 'j':'ⱼ', 'k':'ₖ', 'l':'ₗ', 'm':'ₘ', 'n':'ₙ',
    'o':'ₒ', 'p':'ₚ', 'q':'૧', 'r':'ᵣ', 's':'ₛ', 't':'ₜ', 'u':'ᵤ',
    'v':'ᵥ', 'w':'w', 'x':'ₓ', 'y':'ᵧ', 'z':'₂',

    '0':'₀', '1':'₁', '2':'₂', '3':'₃', '4':'₄',
    '5':'₅', '6':'₆', '7':'₇', '8':'₈', '9':'₉',

    '+':'₊', '-':'₋', '=':'₌', '(':'₍', ')':'₎'
})

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

def align_xy(X, y):
    X = np.asarray(X)
    y = np.asarray(y).ravel()
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    mask = (np.isfinite(y) & np.isfinite(X).all(axis=1))
    X = X[mask]
    y = y[mask]
    if len(y) == 0:
        raise ValueError('No finite samples remain after alignment.')
    return X, y

# --------------- text helpers --------------- #
def rename(name: str | list[str], sep: str | None = '_',
           prefix: str | None = None, suffix: str | None = None):
    
    sep = '' if sep is None else sep
    def _rename(x):
        if suffix is not None:
            x += f'{sep}{suffix}'
        if prefix is not None:
            x = f'{prefix}{sep}{x}'
        return x
    
    if isinstance(name, str):
        return _rename(name)
    return [_rename(x) for x in name]

def lett_encode(n: int):
    label = ''
    while True:
        n, rem = divmod(n, 26)
        label = chr(65 + rem) + label
        if n == 0:
            break
        n -= 1
    return label
    
def superscript(text: str):
    return str(text).translate(_SUPER_MAP)

def subscript(text: str):
    return str(text).translate(_SUB_MAP)

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


# --------------- generators --------------- #
def iter_shermo(lines):
    block = []
    for line in lines:
        if line.startswith('====='):
            if block:
                yield block
                block = []
        block.append(line.strip())
    if block:
        yield block


# --------------- refactor latter --------------- #
"""
import os, shutil, glob
from decimal import Decimal

################################################################################ Text editing

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



################################################################################ OTHERS
def energy_from_xyz(path):
    with open(path) as file:
        comment = file.readlines()[1]
        energy = re.findall(float_regex, comment)[0]
        return float(energy)
"""
