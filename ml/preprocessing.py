from dataclasses import dataclass, asdict
from pathlib import Path

import json
import joblib
import numpy as np
import pandas as pd

from scipy.linalg import qr

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


# ----- preprocessing main config -----
@dataclass(slots=True)
class PreprocessConfig:
    drop_na: bool | None = True
    drop_zero_varr: bool | None = True
    drop_dependent: bool | None = True
    drop_corr: bool | None = True
    scale_data: bool | None = True
    pca: bool | None = False
    # Settings
    linear_tolerance: float = 1e-8
    corr_threshold: float = 0.95
    transformer: object | None = StandardScaler()
    pca_components: int | float | None = 0.95
    pca_solver: str = 'full'


# ----- main transformers -----
class MissingValuesRemover(BaseEstimator, TransformerMixin):
    def __init__(self):
        pass
        
    def fit(self, X, y=None):
        self.remove_columns_ = X.columns[X.isna().sum() > 0]
        return self
        
    def transform(self, X):
        return X.drop(columns=self.remove_columns_, errors='ignore')

class ZeroVarianceRemover(BaseEstimator, TransformerMixin):
    def __init__(self):
        pass
        
    def fit(self, X, y=None):
        self.remove_columns_ = X.columns[X.std() == 0]
        return self
        
    def transform(self, X):
        return X.drop(columns=self.remove_columns_, errors='ignore')

class LinearDependencyRemover(BaseEstimator, TransformerMixin):
    def __init__(self, tolerance: float = 1e-8):
        self.tolerance = tolerance
        
    def fit(self, X, y=None):
        Q, R, P = qr(X, pivoting=True)
        diag = np.abs(np.diag(R))
        n_dependent = np.sum(diag < self.tolerance)
        if n_dependent == 0:
            self.remove_columns_ = []
        else:
            self.remove_columns_ = (X.columns[P[-n_dependent:]].tolist())
        return self
    
    def transform(self, X):
        return X.drop(columns=self.remove_columns_, errors='ignore')

class CorrelationRemover(BaseEstimator, TransformerMixin):
    def __init__(self, threshold: float = 0.95):
        self.threshold = threshold
        
    def fit(self, X, y=None):
        corr_matrix = X.corr().abs()
        upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
        self.remove_columns_ = [column for column in upper.columns if any(upper[column] > self.threshold)]
        return self

    def transform(self, X):
        return X.drop(columns=self.remove_columns_, errors='ignore')

class DataScaler(BaseEstimator, TransformerMixin):
    def __init__(self, scaler=None):
        self.scaler = scaler or StandardScaler()
        
    def fit(self, X, y=None):
        self.scaler.fit(X)
        self.columns_ = X.columns
        self.index_name_ = X.index.name
        return self

    def transform(self, X):
        X_scaled = self.scaler.transform(X)
        return pd.DataFrame(X_scaled, index=X.index, columns=self.columns_)

class PCATransformer(BaseEstimator, TransformerMixin):
    def __init__(self, n_components=0.95, pca_solver='full'):
        self.n_components = n_components
        self.pca_solver = pca_solver

    def fit(self, X, y=None):
        self.pca_ = PCA(n_components=self.n_components, svd_solver=self.pca_solver)
        self.pca_.fit(X)
        return self

    def transform(self, X):
        X_pca = self.pca_.transform(X)
        columns = [f'PC{i + 1}' for i in range(X_pca.shape[1])]
        return pd.DataFrame(X_pca, index=X.index, columns=columns)


# ----- main helpers -----
def save_config(config: PreprocessConfig, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    config_dict = asdict(config)
    if config.transformer is not None:
        config_dict['transformer'] = (config.transformer.__class__.__name__)
    with open(path, 'w') as f:
        json.dump(config_dict, f, indent=4)

def build_pipeline(config: PreprocessConfig) -> Pipeline:
    steps = []
    if config.drop_na:
        steps.append(('remove_missing_values', MissingValuesRemover()))
    if config.drop_zero_varr:
        steps.append(('remove_zero_variance', ZeroVarianceRemover()))
    if config.drop_dependent:
        steps.append(('linear_dependent_columns', LinearDependencyRemover(tolerance=config.linear_tolerance)))
    if config.drop_corr:
        steps.append(('highlly_correlated_columns', CorrelationRemover(threshold=config.corr_threshold)))
    if config.scale_data:
        steps.append(('scale_transform', DataScaler(scaler=config.transformer)))
    if config.pca:
        steps.append(('pca_transform', PCATransformer(n_components=config.pca_components, pca_solver=config.pca_solver)))
    return Pipeline(steps)

def save_pipeline(pipeline: Pipeline, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(pipeline, path)

def load_pipeline(path: str | Path) -> Pipeline:
    return joblib.load(path)

# ----- API ----
def preprocess_data(data: pd.DataFrame, config: PreprocessConfig):
    pipeline = build_pipeline(config)
    transformed = pipeline.fit_transform(data)
    return transformed, pipeline