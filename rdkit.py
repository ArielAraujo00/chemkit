# Import libraries
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from pathlib import Path

# Silence RDKit (AMÉM DEUS)
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from IPython.display import SVG, display
from openbabel import openbabel
from tqdm.auto import tqdm

### UP TO DATE
def convert_MolToXYZ(mol, path='.', reorder=False, mapper=None):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    out = path / f'{mol.GetProp("_Name")}.xyz'
    with out.open('w') as xyz:
        xyz_block = Chem.MolToXYZBlock(mol)
        if not reorder:
            xyz.write(xyz_block)
        else:
            if mapper is not None:
                lines = xyz_block.split('\n')
                header = lines[:2]
                reorder = [lines[v+2] for k, v in mapper.items()]
                rest = [l for l in lines[2:] if not l in reorder]
                xyz_block = "\n".join(header + reorder + rest)
                xyz.writelines(xyz_block)

### NEEDS REVIEW
def read_sdf(file_path):
    return [mol for mol in Chem.SDMolSupplier(file_path)]

def read_txt(file_path):
    with open(file_path) as handle:
        return [Chem.MolFromSmiles(smiles.strip()) for smiles in handle]

def create_filter(mols_series, smarts_patt):
    patt = Chem.MolFromSmarts(smarts_patt)
    return [mol.HasSubstructMatch(patt) for mol in mols_series]

def merge_filters(filters_df):
    return [row[1].any() for row in filters_df.iterrows()]

def get_matches(mols_series, smarts_patt):
    patt = Chem.MolFromSmarts(smarts_patt)
    return [mol.GetSubstructMatches(patt) for mol in mols_series]

def display_molecules(rdmol_series, legend_series=None, prt=None,
                      highlight_query=None, highlight_atoms=None, highlight_bonds=None,
                      n_cols=5, mol_height=340, mol_width=240, Hs=False):
    
    assert isinstance(rdmol_series, pd.core.series.Series) or rdmol_series is None, 'rdmol_series should be a pd.Series'
    assert isinstance(legend_series, pd.core.series.Series) or legend_series is None, 'legend_series should be a pd.Series'
    assert isinstance(prt, str) or prt is None, 'prt should be a path string'
    assert isinstance(highlight_query, str) or highlight_query is None, 'highlight_query shoud be a SMARTS string'
    assert isinstance(highlight_atoms, pd.core.series.Series) or highlight_atoms is None, 'highlight_query should be a pd.Series'
    assert isinstance(highlight_bonds, list) or highlight_bonds is None, 'highlight_bonds should be a list of atom_index'
    
    # Page config
    n_mols = len(rdmol_series.index)
    n_lines = int(np.ceil(n_mols / n_cols))
    page_height = mol_height * n_lines
    page_width = mol_width * n_cols

    # Display config
    d2d = rdMolDraw2D.MolDraw2DSVG(page_width, page_height, mol_width, mol_height)
    opts = d2d.drawOptions()
    opts.legendFraction = 0.25
    opts.legendFontSize = 16
    opts.minFontSize = 14
    opts.drawMolsSameScale = False
    opts.fixedBondLength = 20
    opts.prepareMolsBeforeDrawing = True

    # Highlight (From Query)
    if highlight_query is not None:
        highlight_query = Chem.MolFromSmarts(highlight_query)
        highlight_atom_lists = list()
        highlight_bond_lists = list()
        for mol in rdmol_series:
            atom_matches = mol.GetSubstructMatches(highlight_query)
            atoms_to_highlight = set()
            bonds_to_highlight = set()
            for match in atom_matches:
                atoms_to_highlight.update(match)
                for bond in highlight_query.GetBonds():
                    b = mol.GetBondBetweenAtoms(match[bond.GetBeginAtomIdx()], match[bond.GetEndAtomIdx()])
                    if b is not None:
                        bonds_to_highlight.add(b.GetIdx())
            highlight_atom_lists.append(list(atoms_to_highlight))
            highlight_bond_lists.append(list(bonds_to_highlight))
    else:
        highlight_atom_lists = [[] for _ in rdmol_series]
        highlight_bond_lists = [[] for _ in rdmol_series]

    # Highligh Atoms
    if highlight_atoms is not None:
        highlight_atom_lists = list()
        for atoms in highlight_atoms:
            atoms_to_highlight = set()
            for match in atoms:
                atoms_to_highlight.update(match)
            highlight_atom_lists.append(list(atoms_to_highlight))
    
    # Highligh Bonds
    if highlight_bonds is not None:
        highlight_bond_lists = list()
        for mol, atoms in zip(rdmol_series, highlight_atoms):
            bonds_to_highlight = set()
            for match in atoms:
                for start, end in highlight_bonds:
                    b = mol.GetBondBetweenAtoms(match[start], match[end])
                    if b is not None:
                        bonds_to_highlight.add(b.GetIdx())
            highlight_bond_lists.append(list(bonds_to_highlight))
    
    # Hydrogens
    if Hs:
        rdmol_series = rdmol_series.apply(Chem.AddHs)
    else:
        rdmol_series = rdmol_series.apply(Chem.RemoveHs)

    # Legends
    if legend_series is not None:
        legend_series = list(legend_series)
    else:
        legend_series = [str(i) for i in rdmol_series.index]
    
    # Draw molecules with highlighting and legends
    d2d.DrawMolecules(list(rdmol_series),
                      highlight_atom_lists,
                      highlight_bond_lists,
                      None,  # highlightAtomColors
                      None,  # highlightBondColors
                      None,  # highlightAtomRadii
                      None,  # confIds
                      legend_series)
    
    d2d.FinishDrawing()
    svg = d2d.GetDrawingText()

    # Output
    if prt is not None:
        with open(prt, mode='w') as file:
            file.write(svg)
        return
    return display(SVG(svg))

def set_name(mol_series, legend_series):
    for i in tqdm(mol_series.index):
        mol_series[i].SetProp('_Name', legend_series[i])
    
def embed_molecule(mol, maxiter=1000, prt=False):
    # Default Embedder -> Experimental-Torsion Knowledge Distance Geometry v3 (most recent)
    params = Chem.rdDistGeom.ETKDGv3()
    params.maxIterations=maxiter
    params.randomSeed=42
    success = AllChem.EmbedMolecule(mol, params)

    # If Failure try new aproachs
    if success != 0:
        if prt:
            print('Embedding Failure...\nUpdating params -> useBasicKnowledge=True, useSmallRingTorsions=True, useMacrocycleTorsions=True')
        params.useBasicKnowledge = True
        params.useSmallRingTorsions = True
        params.useMacrocycleTorsions = True
        success = AllChem.EmbedMolecule(mol, params)
    
    if success != 0:
        if prt:
            print('Embedding Failure...\nUpdating params -> ignoreSmoothingFailures=True')
        params.ignoreSmoothingFailures = True
        success = AllChem.EmbedMolecule(mol, params)
    
    if success != 0:
        if prt:
            print('Embedding Failure...\nUpdating params -> enforceChirality=False')
        params.enforceChirality = False
        success = AllChem.EmbedMolecule(mol, params)
    
    if success != 0:
        if prt:
            print('Embedding Failure...\nUpdating params -> useRandomCoords=True')
        params.useRandomCoords = True
        success = AllChem.EmbedMolecule(mol, params)

    if prt:
        print('Sccessfull Embeding...')

def convert_3d(mol_series):
    # To convert molecule to 3d geometry:
    # Add Hydrogens -> Compute 2d Coord -> Compute 3d Coord (Embed) -> Optmize geometry (Optional)
    mol_series = mol_series.apply(Chem.AddHs)
    for i in tqdm(mol_series.index):
        AllChem.Compute2DCoords(mol_series[i])
        embed_molecule(mol_series[i])
    return mol_series

def optMMFF(mol_series, ff_variant='MMFF94', maxiter=1000, core_df=None, constraints=None):

    if constraints is None:
        constraints = {}

    for i in tqdm(mol_series.index):

        mol = mol_series[i]

        # --- build MMFF ---
        mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=ff_variant)
        if mp is None:
            continue

        ff_obj = AllChem.MMFFGetMoleculeForceField(mol, mp)
        if ff_obj is None:
            continue

        conf = mol.GetConformer(0)

        if core_df is not None and constraints:
            row = core_df.loc[i]

            # --- geometry initialization (before minimization) ---
            for a_lbl, b_lbl, d in constraints.get("bond", []):
                a = int(row[a_lbl])
                b = int(row[b_lbl])
                Chem.rdMolTransforms.SetBondLength(conf, a, b, float(d))

            for a_lbl, b_lbl, c_lbl, ang in constraints.get("angle", []):
                a = int(row[a_lbl])
                b = int(row[b_lbl])
                c = int(row[c_lbl])
                Chem.rdMolTransforms.SetAngleDeg(conf, a, b, c, float(ang))

            for a_lbl, b_lbl, c_lbl, d_lbl, dih in constraints.get("dihedral", []):
                a = int(row[a_lbl])
                b = int(row[b_lbl])
                c = int(row[c_lbl])
                d = int(row[d_lbl])
                Chem.rdMolTransforms.SetDihedralDeg(conf, a, b, c, d, float(dih))

            # --- distance constraints (during minimization) ---
            for a_lbl, b_lbl, dist, k in constraints.get("distance", []):
                a = int(row[a_lbl])
                b = int(row[b_lbl])
                ff_obj.AddDistanceConstraint(a, b, float(dist), float(dist), float(k))

        # --- minimize ---
        ff_obj.Minimize(maxIts=maxiter)

def obabel(input_string, input_format, output_format):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(input_format, output_format)
    obabel_mol = openbabel.OBMol()
    obConversion.ReadString(obabel_mol,input_string)
    return obConversion.WriteString(obabel_mol)



def get_atoms_reference(mol_series, smarts, atoms_id):
    smarts = Chem.MolFromSmarts(smarts)
    return pd.DataFrame([mol.GetSubstructMatch(smarts) for mol in mol_series],
                        index=mol_series.index,
                        columns=atoms_id)

def get_prop(pdseries, prop):
    output = list()
    for mol in pdseries:
        try:
            output.append(mol.GetProp(prop))
        except:
            output.append(None)
    return output

def align_structures(mol_series, smiles_str):
    print('\nPlease wait while we align all structures:')
    patt = Chem.MolFromSmarts(smiles_str)
    AllChem.Compute2DCoords(patt)
    for i in tqdm(mol_series.index):
        AllChem.GenerateDepictionMatching2DStructure(mol_series[i], patt)

##################################################################### to implement latter
def conformer_complexity(df, classify=True):
    print('incomplet function...')
    return
    """
    df['sigma_bonds'] = [len(mol.GetSubstructMatches(Chem.MolFromSmarts('[a,A!R]-[!H,!F,!Cl,!Br,!I,$(A#A)]'))) for mol in df['RDKit']]
    df['ring_bonds'] = [len(mol.GetSubstructMatches(Chem.MolFromSmarts('[AR!r3]-[AR!r3]'))) for mol in df['RDKit']]
    if classify:
        df.sort_values(by=['sigma_bonds', 'ring_bonds'], inplace=True)
    """
