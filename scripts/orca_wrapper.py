from pathlib import Path

import pandas as pd

from chemkit.orca import OrcaParser
from chemkit.orca.utils import dump_json

def orca_wrapper(output_path, core_idx, file_id=None, bkp_path=None, overwrite_bkp=False):
    amide_reference = '[#7X3]([#1,#6])([#1,#6])[#6X3](=[#8X1])[#6]'  # Amide SMARTS // Amide Chemspace Project
    smarts_idx = {'N': 0, 'C': 3, 'O': 4, 'R1': 5, 'R2': 1, 'R3': 2}  # Amide atom index from SMARTS != core_idx
    
    # main mapper -> core atoms and its respective labels
    label_map = {core_idx[v]: k for k, v in smarts_idx.items()}

    # Start parser and main variables
    op = OrcaParser(output_path)
    if file_id is None:
        file_id = op.file_path.stem

    # Check for backups
    if bkp_path is None:
        bkp_path = Path.cwd() / f'{file_id}.json'
    else:
        bkp_path = Path(bkp_path) / f'{file_id}.json'
    if bkp_path.is_file() and not overwrite_bkp:
        return

    # Parse descriptors and store results
    op.parse()
    results = {}

    # Molecular Based Descriptors
    if hasattr(op, 'Time'):
        results['Time'] = op.Time
    if hasattr(op, 'SPE'):
        results['EE'] = op.SPE
    if hasattr(op, 'MO'):
        results['HOMO'] = op.MO['occupied'][-1]
        results['LUMO'] = op.MO['virtual'][0]
    if hasattr(op, 'Dipole'):
        results['Dipole'] = op.Dipole['magnitude']
    
    # Atom Based Descriptors
    for idx, atom in label_map.items():
        if hasattr(op, 'Mulliken'):
            results[f'Mulliken_{atom}'] = op.Mulliken[idx]
        if hasattr(op, 'Hirshfeld'):
            results[f'Hirshfeld_{atom}'] = op.Hirshfeld['charge'][idx]
        if hasattr(op, 'Loewdin'):
            results[f'Loewdin_{atom}'] = op.Loewdin[idx]
        if hasattr(op, 'Mayer'):
            results[f'Mayer_Bonded_{atom}'] = op.Mayer['bonded_valence'][idx]
            results[f'Mayer_Free_{atom}'] = op.Mayer['free_valence'][idx]
            results[f'Mayer_Total_{atom}'] = op.Mayer['total_valence'][idx]
        if hasattr(op, 'NPA'):
            results[f'NPA_Charge_{atom}'] = op.NPA['charge'][idx]
            results[f'NPA_Core_{atom}'] = op.NPA['core'][idx]
            results[f'NPA_Valence_{atom}'] = op.NPA['valence'][idx]
            results[f'NPA_Rydberg_{atom}'] = op.NPA['rydberg'][idx]
            results[f'NPA_Total_{atom}'] = op.NPA['total'][idx]
    
    # NBO - Atom and Bond Based Descriptors
    nbo_data = op.get_nbo(atom1=list(core_idx), typ=['LP', 'LV', 'BD', 'BD*'])
    if nbo_data is not None:
        for row in nbo_data.itertuples(index=False):
            l1 = label_map.get(row.atom1)
            if pd.isna(row.atom2):
                key = f'NBO_{row.type}{row.order}_{l1}'
            else:
                l2 = label_map.get(row.atom2)
                key = f'NBO_{row.type}{row.order}_{l1}_{l2}'
            results[key + '_energy'] = row.energy
            results[key + '_occ'] = row.occ
    
    # Save and exit
    dump_json(bkp_path, results)