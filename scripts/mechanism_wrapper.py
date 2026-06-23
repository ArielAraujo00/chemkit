import pandas as pd
from pathlib import Path

from chemkit.orca import OrcaParser
from chemkit import core

def parse_name_mechanism(file):
    name = file.stem
    group = file.parent.parent.stem
    code, soft, calc, lvth = name.rsplit('_', 3)
    conf = 0
    if '_' in code:
        code, conf = code.split('_')
        conf = int(conf)
    return {'file': name, 'group': group, 'code': code, 'conf': conf,
            'soft': soft, 'calc': calc, 'lvth': lvth, 'path': file}

def parse_name_thermo(file):
    name = file.stem
    code, soft, calc, lvth = name.rsplit('_', 3)
    group, mol_id = code.split('-')
    return {'file': name, 'group': group, 'code': code, 'BA': f'BA-{mol_id}',
            'soft': soft, 'calc': calc, 'lvth': lvth, 'path': file}
    
def mechanism_wrapper(files, name_parser):
    data = []
    for f in files:
        op = OrcaParser(f)
        op.parse()
        descript = name_parser(f)
        descript['freq'] = op.Vibrational['imaginary'][0] if op.Vibrational['imaginary'] else None
        descript['energy'] = op.SPE
        descript['gibbs'] = op.Energies['energies']['G']
        data.append(descript)
    return pd.DataFrame(data)