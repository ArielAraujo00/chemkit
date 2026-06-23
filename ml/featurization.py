import numpy as np
import pandas as pd

def get_conformers(index):
    split = index.to_series().str.rsplit('_', n=1, expand=True)
    if split.shape[1] != 2:
        raise ValueError('Index entries must contain at least one "_" separator.')
    split.columns = ['molecule_name', 'conformer_id']
    return split

def boltzmann_weights(series, R=0.008314462618, T=298.15):
    conformers = get_conformers(series.index)
    grouped = series.groupby(conformers['molecule_name'])
    delta_e = grouped.transform(lambda x: x - x.min())
    weights = np.exp(-delta_e / (R * T))
    return weights / grouped.transform(lambda x: np.exp(-(x - x.min()) / (R * T)).sum())

def _rename_columns(df, prefix=None):
    if prefix is None:
        return df
    out = df.copy()
    out.columns = [f'{prefix}_{col}' for col in out.columns]
    return out


def dist_boltzmann(df, reference_column, prefix=None):
    conformers = get_conformers(df.index)
    weights = boltzmann_weights(df[reference_column])
    weighted = df.drop(columns=reference_column).multiply(weights, axis=0)
    dist = weighted.groupby(conformers['molecule_name']).sum()
    dist.index.name = None
    return _rename_columns(dist, prefix)


def dist_select_conformer(df, reference_column, mode='min', prefix=None):
    conformers = get_conformers(df.index)
    grouped = df.groupby(conformers['molecule_name'])[reference_column]
    reference_values = grouped.transform(mode)
    dist = df.drop(columns=reference_column).loc[df[reference_column] == reference_values].copy()
    dist.index = dist.index.str.rsplit('_', n=1).str[0]
    dist.index.name = None
    return _rename_columns(dist, prefix)


def dist_statistical(df, mode='mean', prefix=None):
    conformers = get_conformers(df.index)
    dist = (df.groupby(conformers['molecule_name']).agg(mode))
    dist.index.name = None
    return _rename_columns(dist, prefix)