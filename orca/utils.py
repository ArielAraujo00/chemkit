from tqdm import tqdm
import pandas as pd

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

# --------------- patterns --------------- #
# _FLOAT_RE = re.compile(r'[-]?[0-9]+\.+[0-9]+') # Old version I was using...
_INT_RE   = re.compile(r'[+-]?\d+')
_FLOAT_RE = re.compile(r'[+-]?\d+\.\d+')
_SCY_RE   = re.compile(r'[+-]?\d*\.\d+[Ee][+-]?\d+')
_NUMBER_RE = re.compile(f'{_SCY_RE.pattern}|{_FLOAT_RE.pattern}|{_INT_RE.pattern}')

_NBO_RE = re.compile(
    r"""
    ^\s*\d+\.\s+                             # index
    (?P<type>CR|LP|LV|BD\*?|RY)\s*           # NBO type
    \(\s*(?P<order>\d+)\s*\)\s+
    (?:
        (?P<elem1>[A-Za-z]+)\s*(?P<a1>\d+)   # first atom
        (?:\s*-\s*
            (?P<elem2>[A-Za-z]+)\s*(?P<a2>\d+)
        )?
    )
    \s+
    (?P<occ>[-+]?\d*\.\d+)\s+                # occupancy
    (?P<E>[-+]?\d*\.\d+)                     # energy
    """,
    re.VERBOSE,
)

# --------------- helpers --------------- #
def _parse_path(file_path):
    assert isinstance(file_path, (str, Path)), 'Must be a path-like object (str | Path)'
    if not isinstance(file_path, Path):
        return Path(file_path)
    else:
        return file_path

def _read_lines(file):
    file = _parse_path(file)
    with file.open('r') as f:
        return f.readlines()

def _extract_numbers(line):
    out = []
    for s in _NUMBER_RE.findall(line):
        if 'e' in s or 'E' in s:
            out.append(float(s))
        elif '.' in s:
            out.append(float(s))
        else:
            out.append(int(s))
    return out

def _safe_float(number):
    try:
        return float(number)
    except (ValueError, TypeError):
        return None

def _safe_int(number):
    try:
        return int(number)
    except (ValueError, TypeError):
        return None

def _seek_tag(lines, tags, start=0):
    assert isinstance(tags, tuple), 'tags must be a tuple (did you use @block?)'
    for i in range(start, len(lines)):
        line = lines[i]
        for tag in tags:
            if tag in line:
                return i, tag
    return None, None

# --------------- endpoints --------------- #
def _blank(line):
    return not line.strip()

def _not_digit(line):
    s = line.strip()
    return (not s) or (not s[0].isdigit())

def _separator(line):
    s = line.strip()
    return s and s[0] in '-=*'

def _contains(text):
    def _fn(line):
        return text in line
    return _fn

# --------------- decorators --------------- #
def block(tags, key, header=0, endpoint=_blank, multiple=False, mode='last'):
    
    # normalize tag -> tuple[str]
    if isinstance(tags, str):
        tags = (tags,)
    else:
        tags = tuple(tags)
    
    # validate mode
    if mode not in {'first', 'last', 'all'}:
        raise ValueError(f"Invalid mode '{mode}'")

    def wrapper(func):
        func.tags = tags
        func.key = key
        func.header = header
        func.endpoint = endpoint
        func.multiple = multiple
        func.mode = mode
        return func
    return wrapper
 
# --------------- parsers --------------- #
def parse_energies(energies, unit='kcal/mol'):
    parsed = {'energies': {'EE': energies['Electronic energy'],
                           'ZPE': (energies['Electronic energy'] + energies['Zero point energy']
                                   if energies['Electronic energy'] is not None
                                   and energies['Zero point energy'] is not None
                                   else None),
                           'TE': energies['Total thermal energy'],
                           'H': energies['Total Enthalpy'],
                           'S': energies['Final entropy term'],
                           'G': energies['Final Gibbs free energy']},
              
              'corrections': {'Zero point correction': energies['Zero point energy'],
                              'Vibrational correction': energies['Thermal vibrational correction'],
                              'Rotational correction': energies['Thermal rotational correction'],
                              'Translational correction': energies['Thermal translational correction'],
                              'Enthalpy correction': energies['Thermal Enthalpy correction'],
                              'Electronic entropy': energies['Electronic entropy'],
                              'Vibrational entropy': energies['Vibrational entropy'],
                              'Rotational entropy': energies['Rotational entropy'],
                              'Translational entropy': energies['Translational entropy']}
             }
    
    for p in parsed:
        parsed[p] = {k: (v * _EH_TO[unit] if v is not None else None)
                     for k, v in parsed[p].items()}

    return parsed

def parse_ir(ir_spectra, eps=1e-3):
    modes = {'real': [], 'imaginary': []}
    for ir in ir_spectra:
        if ir > eps:
            modes['real'].append(ir)
        elif ir < -eps:
            modes['imaginary'].append(ir)
    return modes

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
