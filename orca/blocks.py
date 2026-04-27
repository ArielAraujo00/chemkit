# Built-in
import re

# This Project
from . import utils

# ----- orca parser blocks -----
@block(tags='FINAL SINGLE POINT ENERGY', key='SPE', multiple=True)
def spe(lines, i=0):
    output = None
    i, _ = _seek_tag(lines, spe.tags, start=i)
    if i is None:
        return output
    i += spe.header
    vals = extract_numbers(lines[i])
    output.append(vals[-1])
    return output

@block(tags='ORBITAL ENERGIES', key='MO', header=4, endpoint=_not_digit)
def orbitals(lines, i=0):
    output = {'occupied': [], 'virtual': []}
    i, _ = _seek_tag(lines, orbitals.tags, start=i)
    if i is None:
        return output
    i += orbitals.header
    while i < len(lines) and not orbitals.endpoint(lines[i]):
        vals = extract_numbers(lines[i])
        if vals[1] > 0:
            output['occupied'].append(vals[2])
        else:
            output['virtual'].append(vals[2])
        i += 1
    return output

@block(tags='MULLIKEN ATOMIC CHARGES', key='Mulliken', header=2,
       endpoint=_contains('Sum of atomic charges'))
def mulliken(lines, i=0):
    output = []
    i, _ = _seek_tag(lines, mulliken.tags, start=i)
    if i is None:
        return output
    i += mulliken.header
    while i < len(lines) and not mulliken.endpoint(lines[i]):
        vals = extract_numbers(lines[i])
        output.append(vals[-1])
        i += 1
    return output

@block(tags='LOEWDIN ATOMIC CHARGES', key='Loewdin', header=2, endpoint=_blank)
def loewdin(lines, i=0):
    output = []
    i, _ = _seek_tag(lines, loewdin.tags, start=i)
    if i is None:
        return output
    i += loewdin.header
    while i < len(lines) and not loewdin.endpoint(lines[i]):
        vals = extract_numbers(lines[i])
        output.append(vals[-1])
        i += 1
    return output

@block(tags='MAYER POPULATION ANALYSIS', key='Mayer', header=11, endpoint=_blank)
def mayer(lines, i=0):
    output = {'total_valence': [], 'bonded_valence': [], 'free_valence': []}
    i, _ = _seek_tag(lines, mayer.tags, start=i)
    if i is None:
        return output
    i += mayer.header
    while i < len(lines) and not mayer.endpoint(lines[i]):
        vals = extract_numbers(lines[i])
        output['total_valence'].append(vals[-3])
        output['bonded_valence'].append(vals[-2])
        output['free_valence'].append(vals[-1])
        i += 1
    return output

@block(tags='HIRSHFELD ANALYSIS', key='Hirshfeld', header=7, endpoint=_blank)
def hirshfeld(lines, i=0):
    output = {'charge': [], 'spin': []}
    i, _ = _seek_tag(lines, hirshfeld.tags, start=i)
    if i is None:
        return output
    i += hirshfeld.header
    while i < len(lines) and not hirshfeld.endpoint(lines[i]):
        vals = extract_numbers(lines[i])
        if len(vals) >= 2:
            output['charge'].append(vals[-2])
            output['spin'].append(vals[-1])
        i += 1
    return output

@block(tags='Summary of Natural Population Analysis:', key='NPA',
       header=6, endpoint=_separator)
def npa(lines, i=0):
    output = {'charge': [], 'core': [], 'valence': [], 'rydberg': [], 'total': []}
    i, _ = _seek_tag(lines, npa.tags, start=i)
    if i is None:
        return output
    i += npa.header
    while i < len(lines) and not npa.endpoint(lines[i]):
        vals = extract_numbers(lines[i])
        output['charge'].append(vals[-5])
        output['core'].append(vals[-4])
        output['valence'].append(vals[-3])
        output['rydberg'].append(vals[-2])
        output['total'].append(vals[-1])
        i += 1
    return output

@block(tags='NATURAL BOND ORBITALS (Summary):', key='NBO', header=1,
       endpoint=_contains('NBO analysis completed'))
def nbo(lines, i=0):
    output = {'atom': {'CR': {}, 'LP': {}, 'LV': {}, 'RY': {}},
              'bond': {'BD': {}, 'BD*': {}}}
    i, _ = _seek_tag(lines, nbo.tags, start=i)
    if i is None:
        return output
    i += nbo.header
    while i < len(lines) and not nbo.endpoint(lines[i]):
        line = lines[i]
        m = _NBO_RE.match(line)
        if m:
            typ = m.group('type')
            order = int(m.group('order'))
            occ = _safe_float(m.group('occ'))
            ene = _safe_float(m.group('E'))
            a1 = int(m.group('a1'))
            a2 = m.group('a2')
            if a2 is None:
                key = (a1,)
                target = 'atom'
            else:
                key = (a1, int(a2))
                target = 'bond'
            d = output[target][typ].setdefault(order, {})
            d[key] = {'occupancy': occ, 'energy': ene}
        i += 1
    return output

@block(tags='DIPOLE MOMENT', key='Dipole', endpoint=_contains('Magnitude (Debye)'))
def dipole(lines, i=0):
    output = {'vector': None, 'magnitude': None}
    i, _ = _seek_tag(lines, dipole.tags, start=i)
    if i is None:
        return output
    i += dipole.header
    while i < len(lines) and not dipole.endpoint(lines[i]):
        if 'Total Dipole Moment' in lines[i]:
            vals = extract_numbers(lines[i])
            if len(vals) >= 3:
                output['vector'] = tuple(vals[:3])
        i += 1
    # final line contains magnitude
    vals = extract_numbers(lines[i])
    if vals:
        output['magnitude'] = vals[0]
    return output

@block(tags='Rotational spectrum', key='Rotational', endpoint=_contains('Rotational constants in MHz'))
def rotational_constants(lines, i=0):
    output = {'cm-1': None, 'MHz': None}
    i, _ = _seek_tag(lines, rotational_constants.tags, start=i)
    if i is None:
        return output
    i += rotational_constants.header
    while i < len(lines) and not rotational_constants.endpoint(lines[i]):
        if 'Rotational constants in cm-1' in lines[i]:
            vals = extract_numbers(lines[i])
            if vals:
                output['cm-1'] = tuple(vals)
        i += 1
    # final line contains (MHz)
    vals = extract_numbers(lines[i])
    if vals:
        output['MHz'] = tuple(vals)
    return output

@block(tags='VIBRATIONAL FREQUENCIES', key='Vibrational', header=5, endpoint=_blank)
def ir_spectra(lines, i=0):
    output = []
    i, _ = _seek_tag(lines, ir_spectra.tags, start=i)
    if i is None:
        return output
    i += ir_spectra.header
    while i < len(lines) and not ir_spectra.endpoint(lines[i]):
        vals = extract_numbers(lines[i])
        output.append(vals[0])
        i += 1
    return parse_ir(output)

@block(tags='THERMOCHEMISTRY AT', key='Energies')
def energies(lines, i=0, conv=True, const='kcal/mol'):
    output = {
        'Electronic energy': None,
        'Total thermal energy': None,
        'Total Enthalpy': None,
        'Final entropy term': None,
        'Final Gibbs free energy': None,

        'Zero point energy': None,
        'Thermal vibrational correction': None,
        'Thermal rotational correction': None,
        'Thermal translational correction': None,
        'Thermal Enthalpy correction': None,
        'Electronic entropy': None,
        'Vibrational entropy': None,
        'Rotational entropy': None,
        'Translational entropy': None
    }
    i, _ = _seek_tag(lines, energies.tags, start=i)
    if i is None:
        return output
    i += energies.header
    while i < len(lines):
        for k in output:
            if k in lines[i]:
                vals = extract_numbers(lines[i])
                if vals:
                    output[k] = vals[0]
        i += 1
        if all(v is not None for v in output.values()):
            return parse_energies(output, conv=conv)
    return parse_energies(output, conv=conv, const=const)

