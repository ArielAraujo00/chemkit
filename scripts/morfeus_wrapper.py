from pathlib import Path

import numpy as np
from morfeus import read_xyz, Sterimol, BuriedVolume, Pyramidalization, Dispersion, SASA

from chemkit.orca.utils import dump_json, _safe_float

##################################  MORFEUS  #######################################
# Sterimol
def get_sterimol(elements, coords, atom1, atom2, rtype='crc', encoding=None, py_index=True):
    if encoding is None:
        encoding = f'{elements[atom1]}{atom1}_{elements[atom2]}{atom2}'
    if py_index:
        atom1, atom2 = atom1+1, atom2+1
        
    # Compute Params
    stmol = Sterimol(elements, coords, atom1, atom2, radii_type=rtype)

    # Parse Output
    output = {f'Sterimol_B1_{encoding}': stmol.B_1_value,
              f'Sterimol_B5_{encoding}': stmol.B_5_value,
              f'Sterimol_L_{encoding}': stmol.L_value}
    
    return output

# SterimolBur
def get_sterimol_bur(elements, coords, atom1, atom2, rtype='crc', bur_method='delete',
                     r=5.5, scale=0.5, encoding=None, py_index=True):
    if encoding is None:
        encoding = f'{elements[atom1]}{atom1}_{elements[atom2]}{atom2}'
    if py_index:
        atom1, atom2 = atom1+1, atom2+1

    # Compute Params
    stmol = Sterimol(elements, coords, atom1, atom2, radii_type=rtype)
    stmol.bury(method=bur_method, sphere_radius=r, radii_scale=scale)

    # Parse Output
    output = {f'SterimolBur_{r:.1f}A_B1_{encoding}': stmol.B_1_value,
              f'SterimolBur_{r:.1f}A_B5_{encoding}': stmol.B_5_value,
              f'SterimolBur_{r:.1f}A_L_{encoding}': stmol.L_value}
    
    return output

# VBur
def get_vbur(elements, coords, atom1, r, rtype='bondi', rscale=1.17,
             Hs=False, encoding=None, py_index=True):
    if encoding is None:
        encoding = f'{elements[atom1]}{atom1}'
    if py_index:
        atom1 = atom1+1
    
    # Compute Params
    bv = BuriedVolume(elements, coords, atom1, radius=r, include_hs=Hs,
                      radii_type=rtype, radii_scale=rscale)

    # Parse Output
    output = {f'VBur_{r:.1f}A_bury_{encoding}': bv.buried_volume,
              f'VBur_{r:.1f}A_free_{encoding}': bv.free_volume,
              f'VBur_{r:.1f}A_frac_{encoding}': bv.fraction_buried_volume}
    
    return output

# VBur Sectors
def get_vbur_sectors(elements, coords, atom1, r, zaxis, xzplane, rtype='bondi',
                     rscale=1.17, Hs=False, encoding=None, py_index=True):
    if encoding is None:
        encoding = f'{elements[atom1]}{atom1}'
    if py_index:
        atom1 = atom1+1
    
    # Compute Params
    bv = BuriedVolume(elements, coords, atom1, radius=r, include_hs=Hs, radii_type=rtype,
                      radii_scale=rscale, z_axis_atoms=zaxis, xz_plane_atoms=xzplane)
    bv.octant_analysis()

    # Parse Output
    output = {}
    space_maps = {'qrd': bv.quadrants,
                  'oct': bv.octants}
    metric_maps = {'bury': 'buried_volume',
                   'free': 'free_volume',
                   'frac': 'percent_buried_volume'}
    for space_code, space_data in space_maps.items():
        for metric_code, metric_key in metric_maps.items():
            for idx, value in space_data[metric_key].items():
                output[f"VBur_{r:.1f}A_{space_code}{idx}_{metric_code}_{encoding}"] = value
    
    return output

# SASA
def get_sasa(elements, coords, probe=1.4, atoms=None, rtype='crc', encoding=None, py_index=True):
    if atoms is not None:
        atoms = list(atoms)
        
        if encoding is None:
            encoding = [f'{elements[a]}{a}' for a in atoms]
        else:
            encoding = list(encoding)
        
        if py_index:
            atoms = [a+1 for a in atoms]

    # Compute Parameters:
    sasa = SASA(elements, coords, probe_radius=probe, radii_type=rtype)

    # Parse output
    output = {f'SASA_{probe:.1f}A_total_area': sasa.area,
              f'SASA_{probe:.1f}A_total_volume': sasa.volume}
    if atoms is not None:
        output.update({f'SASA_{probe:.1f}A_{label}_area': sasa.atom_areas[idx] for label, idx in zip(encoding, atoms)})

    return output

# Dispersion
def get_dispersion(elements, coords, atoms=None, rtype='rahm', encoding=None, py_index=True):
    if atoms is not None:
        atoms = list(atoms)
        
        if encoding is None:
            encoding = [f'{elements[a]}{a}' for a in atoms]
        else:
            encoding = list(encoding)
        
        if py_index:
            atoms = [a+1 for a in atoms]
    
    # Compute Parameters:
    disp = Dispersion(elements, coords, radii_type=rtype)
    
    # Parse output
    output = {f'Dispersion_total_area': disp.area,
              f'Dispersion_total_volume': disp.volume,
              f'Dispersion_Pint': disp.p_int,
              f'Dispersion_Pmax': disp.p_max,
              f'Dispersion_Pmin': disp.p_min}
    
    if atoms is not None:
        output.update({f'Dispersion_area_{label}': disp.atom_areas[idx] for label, idx in zip(encoding, atoms)})
        output.update({f'Dispersion_Pint_{label}': disp.atom_p_int[idx] for label, idx in zip(encoding, atoms)})
        output.update({f'Dispersion_Pmax_{label}': disp.atom_p_max[idx] for label, idx in zip(encoding, atoms)})
        output.update({f'Dispersion_Pmin_{label}': disp.atom_p_min[idx] for label, idx in zip(encoding, atoms)})

    return output

# Pyramidization
def get_pyramidization(elements, coords, atom1, neighbors=None, rtype='pyykko',
                       rscale=1.2, method='distance', encoding=None, py_index=True):
    if encoding is None:
        encoding = f'{elements[atom1]}{atom1}'
    if py_index:
        atom1 = atom1+1
        if neighbors is not None:
            neighbors = list(neighbors)
            neighbors = [n+1 for n in neighbors]
            
            
    
    # Compute Params
    pyr = Pyramidalization(coords, atom1, neighbor_indices=neighbors, elements=elements,
                           radii_type=rtype, method=method, scale_factor=rscale)

    # Parse Output
    output = {f'Pyramidalization_{encoding}_alpha': pyr.alpha,
              f'Pyramidalization_{encoding}_P': pyr.P,
              f'Pyramidalization_{encoding}_P_angle': pyr.P_angle}
    
    return output

##################################  GEOMETRIC  #######################################
# dont recall why i made that for lol
def vector_array_wrapper(func):
    def wrapper(*args, **kwargs):
        new_args = [np.array(arg) for arg in args]
        return func(*new_args, **kwargs)
    return wrapper

@vector_array_wrapper
def normalize_vector(v):
    return v / np.linalg.norm(v)

@vector_array_wrapper
def line_vector(p1, p2):
    return p2 - p1

@vector_array_wrapper
def plane_normal(p1, p2, p3):
    v1 = line_vector(p1, p2)
    v2 = line_vector(p1, p3)
    return np.cross(v1, v2)

@vector_array_wrapper
def distance(p1, p2):
    return np.linalg.norm(p2 - p1)

@vector_array_wrapper
def angle(p1, p2, p3):
    v1 = line_vector(p1, p2)
    v2 = line_vector(p3, p2)
    dot = np.dot(v1, v2)
    norm_prod = np.linalg.norm(v1) * np.linalg.norm(v2)
    cos_theta = np.clip(dot / norm_prod, -1.0, 1.0)
    return np.degrees(np.arccos(cos_theta))

@vector_array_wrapper
def dihedral(p1, p2, p3, p4):
    b0 = -line_vector(p1, p2)
    b1 = line_vector(p2, p3)
    b2 = line_vector(p3, p4)
    b1_unit = normalize_vector(b1)
    v = np.cross(b0, b1_unit)
    w = np.cross(b2, b1_unit)
    x = np.dot(v, w)
    y = np.dot(np.cross(v, w), b1_unit)
    return np.degrees(np.arctan2(y, x))

##################################  WRAPPER #######################################

def morfeus_wrapper(xyz_path, core_atoms, file_id=None, bkp_path=None, overwrite_bkp=False):
    amide_reference = '[#7X3]([#1,#6])([#1,#6])[#6X3](=[#8X1])[#6]'  # Amide SMARTS // Amide Chemspace Project
    smarts_idx = {'N': 0, 'C': 3, 'O': 4, 'R1': 5, 'R2': 1, 'R3': 2}  # Amide atom index from SMARTS != core_atoms
    core_bonds = [(0, 1), (0, 2), (0, 3), (3, 4), (3, 5)]  # Amide main bonds
    core_angls = [(0, 3, 4), (4, 3, 5), (5, 3, 0), (1, 0, 3), (3, 0, 2), (2, 0, 1)] # Amide main angles
    core_dihed = [(5, 3, 0, 1), (5, 3, 0, 2), (4, 3, 0, 1), (4, 3, 0, 2)] # Amide main dihedrals

    # main mappers -> core atoms and its respective labels
    label_map = {core_atoms[v]: k for k, v in smarts_idx.items()}
    bonds_idx = [(core_atoms[a1], core_atoms[a2]) for a1, a2 in core_bonds]
    bonds_map = {f'{label_map[a1]}_{label_map[a2]}': (a1, a2) for a1, a2 in bonds_idx}
    angls_idx = [(core_atoms[a1], core_atoms[a2], core_atoms[a3]) for a1, a2, a3 in core_angls]
    angls_map = {f'{label_map[a1]}_{label_map[a2]}_{label_map[a3]}': (a1, a2, a3) for a1, a2, a3 in angls_idx}
    dihed_idx = [(core_atoms[a1], core_atoms[a2], core_atoms[a3], core_atoms[a4]) for a1, a2, a3, a4 in core_dihed]
    dihed_map = {f'{label_map[a1]}_{label_map[a2]}_{label_map[a3]}_{label_map[a4]}': (a1, a2, a3, a4) for a1, a2, a3, a4 in dihed_idx}
    

    # Start parser and main variables ###### rename latter
    if file_id is None:
        file_id = Path(xyz_path).stem

    # Check for backups
    if bkp_path is None:
        bkp_path = Path.cwd() / f'{file_id}.json'
    else:
        bkp_path = Path(bkp_path) / f'{file_id}.json'
    if bkp_path.is_file() and not overwrite_bkp:
        return
        
    # Read xyz and extract descriptors
    elements, coords = read_xyz(xyz_path)
    results = {}
    
    ### Sterimol: All Core Bonds
    for bond, (a1, a2) in bonds_map.items():
        results.update(get_sterimol(elements, coords, a1, a2, encoding=bond))
        results.update(get_sterimol_bur(elements, coords, a1, a2, encoding=bond))

    ### Buried Volume: All Core Atoms | r=3.5 and r=5.0
    for idx, atom in label_map.items():
        results.update(get_vbur(elements, coords, idx, 3.5, encoding=atom))
        results.update(get_vbur(elements, coords, idx, 5.0, encoding=atom))

    ### SASA: Global + Core Atoms | probe=1.4 and r=2.5
    results.update(get_sasa(elements, coords, probe=1.4, atoms=label_map.keys(), encoding=label_map.values()))
    results.update(get_sasa(elements, coords, probe=2.5, atoms=label_map.keys(), encoding=label_map.values()))
    
    ### Dispersion: Global + Core Atoms
    results.update(get_dispersion(elements, coords, atoms=label_map.keys(), encoding=label_map.values()))
    
    ### Pyramidization: Only N
    results.update(get_pyramidization(elements, coords, core_atoms[0], encoding='N'))

    ### Geometric Descriptors
    for dist, (a1, a2) in bonds_map.items():
        results[f'Distance_{dist}'] = distance(coords[a1], coords[a2])
    for angl, (a1, a2, a3) in angls_map.items():
        results[f'Angle_{angl}'] = angle(coords[a1], coords[a2], coords[a3])
    for dhed, (a1, a2, a3, a4) in dihed_map.items():
        results[f'Dihedral_{dhed}'] = dihedral(coords[a1], coords[a2], coords[a3], coords[a4])
    
    # Save and exit
    results = {k: _safe_float(v) for k, v in results.items()}
    dump_json(bkp_path, results)