from tqdm.auto import tqdm
import pandas as pd

import re

# This Project
from chemkit import core

# --------------- patterns --------------- #
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
def seek_tag(lines, tags, start=0):
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
def parse_ir(ir_spectra, eps=1e-3):
    modes = {'real': [], 'imaginary': []}
    for ir in ir_spectra:
        if ir > eps:
            modes['real'].append(ir)
        elif ir < -eps:
            modes['imaginary'].append(ir)
    return modes

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
        parsed[p] = {k: (v * core._EH_TO[unit] if v is not None else None)
                     for k, v in parsed[p].items()}

    return parsed


