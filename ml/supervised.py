import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from adjustText import adjust_text

from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import (r2_score, mean_absolute_error, mean_squared_error, 
                             accuracy_score, precision_score, recall_score,
                             f1_score, roc_auc_score)

from chemkit.core import rename, lett_encode, align_xy
from chemkit.ml.preprocessing import preprocess_data


###################### DATA TABLE ####################################
class Table:
    def __init__(self, data, subsets=None, name=None):
        self.frame = pd.DataFrame(data)
        self.name = name
        self.subsets = subsets or {self.name if self.name is not None else 'main':
                                   self.frame.columns.tolist()}
        self.shape = self.frame.shape
        self.columns = self.frame.columns
        self.index = self.frame.index
    
    def add_subset(self, name, columns):
        if isinstance(columns, str): columns = [columns]
        self.subsets[name] = columns
        return self

    def subset(self, name):
        return self.subsets[name]
    
    def __getitem__(self, key):
        if key not in self.subsets:
            raise KeyError(f'Subset "{key}" not found')
        return TableView(self, self.subsets[key])
    
    def select(self, columns):
        return self.__class__(self.frame[columns], subsets=self.subsets, name=self.name)

    def drop(self, columns):
        return self.__class__(self.frame.drop(columns=columns), subsets=self.subsets, name=self.name)

    def copy(self):
        return self.__class__(self.frame.copy(), subsets=self.subsets.copy(), name=self.name)

    def to_frame(self):
        return self.frame.copy()

    def __repr__(self):
        return (f'{self.__class__.__name__}('
                f'shape={self.shape}, '
                f'subsets={list(self.subsets.keys())}'
                f')')

class FeatureTable(Table):
    def __init__(self, data, subsets=None, name=None):
        super().__init__(data=data, subsets=subsets, name=name)
    
    def X(self):
        return self.frame.select_dtypes(include=np.number)

    def B(self):
        return self.frame.select_dtypes(include=bool)

    def preprocess(self, config):
        transformed, pipeline = preprocess_data(self.X(), config)
        return FeatureTable(transformed, subsets=self.subsets,
                            name=self.name), pipeline
    
    def corr(self):
        return self.X().corr()
    
    def regress(self, targets):
        return screen(self.X(), targets, model='regression')

    def classify(self, targets, regression=False):
        return screen(self.X(), targets, model='classification', regression=regression)

    def __getitem__(self, key):
        if key not in self.subsets:
            raise KeyError(f'Subset "{key}" not found')
        return FeatureView(self, self.subsets[key])

    def __repr__(self):
        return (f'FeatureTable('
                f'shape={self.shape}, '
                f'numeric={self.X().shape[1]}, '
                f'boolean={self.B().shape[1]}, '
                f'subsets={list(self.subsets.keys())}'
                f')')

class TargetTable(Table):
    def __init__(self, data, subsets=None, name=None):
        super().__init__(data=data, subsets=subsets, name=name)

    def Y(self):
        return self.frame.select_dtypes(include=np.number)

    def C(self):
        return self.frame.select_dtypes(include=bool)

    # Logarithm Transform
    def log(self, columns=None, prefix='Log10'):
        return self.lx(10, columns=columns, prefix=prefix)

    def ln(self, columns=None, prefix='Ln'):
        return self.lx(np.e, columns=columns, prefix=prefix)

    def lg(self, columns=None, prefix='Log2'):
        return self.lx(2, columns=columns, prefix=prefix)

    def lx(self, base, columns=None, prefix='LogX', suffix=None, sep='_'):
        frame = self.frame.copy()
        columns = columns or self.Y().columns
        new_columns = rename(columns, prefix=prefix, suffix=suffix, sep=sep)
        if isinstance(columns, str): columns = [columns]
        if isinstance(new_columns, str): new_columwns = [new_columns]    
        for old_col, new_col in zip(columns, new_columns):
            values = frame[old_col].astype(float)
            logged = np.full(len(values), np.nan)
            mask = values > 0
            logged[mask] = np.log(values[mask]) / np.log(base)
            frame[new_col] = logged
        subsets = self.subsets.copy()
        subset_name = prefix if suffix is None else f'{prefix}{sep}{suffix}'
        subsets[subset_name] = new_columns
        return TargetTable(frame, subsets=subsets, name=self.name)

    # Bolean Transform
    def threshold(self, threshold: float, columns=None, inverse=False,
                  prefix='Thr', suffix=None, sep='_'):
        frame = self.frame.copy()
        columns = columns or self.Y().columns
        new_columns = rename(columns, prefix=prefix, suffix=suffix, sep=sep)

        if isinstance(columns, str): columns = [columns]
        if isinstance(new_columns, str): new_columns = [new_columns]
        for old_col, new_col in zip(columns, new_columns):
            frame[new_col] = ((frame[old_col] <= threshold) if inverse
                              else (frame[old_col] >= threshold))
        subsets = self.subsets.copy()
        subset_name = prefix if suffix is None else f'{prefix}{sep}{suffix}'
        subsets[subset_name] = new_columns
        return TargetTable(frame, subsets=subsets, name=self.name)

    # Bin Transform
    def bin(self, bins, columns=None, mode='value', labels=None,
            prefix='Bin', suffix=None, sep='_'):
        frame = self.frame.copy()
        columns = columns or self.Y().columns
        
        if isinstance(columns, str): columns = [columns]
        created_columns = []
        for col in columns:
            values = frame[col].values
            if isinstance(bins, int):
                if mode == 'value':
                    edges = np.linspace(values.min(), values.max(), bins + 1)
                elif mode == 'index':
                    order = np.sort(values)
                    idx = np.linspace(0, len(order) - 1, bins + 1).astype(int)
                    edges = order[idx]
                    edges[0] = -np.inf
                    edges[-1] = np.inf
            else:
                bins = np.asarray(bins, dtype=float)
                bins = bins / bins.sum()
                if mode == 'value':
                    quantiles = np.cumsum(bins)[:-1]
                    edges = [-np.inf, *np.quantile(values, quantiles), np.inf]
                elif mode == 'index':
                    order = np.sort(values)
                    idx = (np.cumsum(bins)[:-1] * len(order)).astype(int)
                    edges = [-np.inf, *order[idx], np.inf]
            n_bins = len(edges) - 1
            current_labels = labels or [lett_encode(i) for i in range(n_bins)]
            new_col = rename(col, prefix=prefix, suffix=suffix, sep=sep)
            frame[new_col] = pd.cut(values, bins=edges, include_lowest=True,
                                    labels=current_labels)
            created_columns.append(new_col)
        subsets = self.subsets.copy()
        subset_name = prefix if suffix is None else f'{prefix}{sep}{suffix}'
        subsets[subset_name] = created_columns
        return TargetTable(frame, subsets=subsets, name=self.name)

    def __getitem__(self, key):
        if key not in self.subsets:
            raise KeyError(f'Subset "{key}" not found')
        return TargetView(self, self.subsets[key])

    def __repr__(self):
        return (f'TargetTable('
                f'shape={self.shape}, '
                f'numeric={self.Y().shape[1]}, '
                f'boolean={self.C().shape[1]}, '
                f'subsets={list(self.subsets.keys())}'
                f')')


###################### TABLE VISUALIZATION API ####################################
class TableView:
    def __init__(self, parent, columns):
        self.parent = parent
        self.columns = columns
    def frame(self):
        return self.parent.frame[self.columns]

class FeatureView(TableView):
    def X(self):
        return self.parent.frame[self.columns].select_dtypes(include=np.number)
    def B(self):
        return self.parent.frame[self.columns].select_dtypes(include=bool)

class TargetView(TableView):
    def Y(self):
        return self.parent.frame[self.columns].select_dtypes(include=np.number)
    def C(self):
        return self.parent.frame[self.columns].select_dtypes(include=bool)


###################### MODELING RESULTS & ANALYSIS ####################################
class ModelBest:
    def __init__(self, results, row):
        self.results = results
        self.row = row
    @property
    def feature(self):
        return self.row['feature']
    @property
    def target(self):
        return self.row['target']
    @property
    def model(self):
        return self.results.models[(self.feature, self.target)]
    @property
    def metrics(self):
        return self.row.drop(['feature', 'target'])
    @property
    def X(self):
        return self.results.features[self.feature]
    @property
    def y(self):
        return self.results.targets[self.target]
    @property
    def y_pred(self):
        X = self.X()
        if X.ndim == 1:
            X = X.to_frame()
        return self.model.predict(self.X())
    @property
    def fit(self):
        X, y = align_xy(self.X(), self.y())
        return pd.DataFrame({'X': X.ravel(), 'y': y,
                             'y_pred': self.model.predict(X)})
    def __repr__(self):
        return (f'ModelBest('
                f'feature={self.feature!r}, '
                f'target={self.target!r})')

class ModelResults:
    def __init__(self, metrics, models=None,
                 features=None, targets=None):
        self.metrics = metrics
        self.models = models or {}
        self.features = features
        self.targets = targets
    def sort(self, metric, n=10, ascending=False):
        return self.metrics.sort_values(metric, ascending=ascending).head(n)
    def filter(self, metric, threshold):
        return self.metrics[self.metrics[metric] >= threshold]
    def best(self, metric, ascending=False, iloc=0):
        row = self.metrics.sort_values(metric, ascending=ascending).iloc[iloc]
        return ModelBest(self, row)
    def keys(self):
        return list(self.models.keys())
    def __repr__(self):
        return (f'{self.__class__.__name__}('
                f'n_models={len(self.models)}, '
                f'metrics={list(self.metrics.columns)}'
                f')')

class RegressionResults(ModelResults):
    def predict(self, feature, target, X):
        return self.models[(feature, target)].predict(X)

class ClassificationResults(ModelResults):
    def __init__(self, metrics, models=None,
                 features=None, targets=None,
                 thresholds=None, regressions=None):
        super().__init__(metrics, models=models)
        self.thresholds = thresholds or {}
        self.regressions = regressions or {}

    def predict(self, feature, target, X):
        return self.models[(feature, target)].predict(X)

    def predict_proba(self, feature, target, X):
        return self.models[(feature, target)].predict_proba(X)

    def threshold(self, feature, target):
        return self.thresholds[(feature, target)]

    def regression(self, feature, target):
        return self.regressions[(feature, target)]


###################### MODELING METRICS ####################################
def regression_metrics(X, y_true, y_pred, n_features=1):
    n = len(y_true)
    dof = n - n_features - 1
    residuals = y_true - y_pred
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    metrics =  {'n': n, 'dof': dof, 'r2': r2_score(y_true, y_pred),
                'mae': mean_absolute_error(y_true, y_pred),
                'rmse': np.sqrt(mean_squared_error(y_true, y_pred)),
                'rss': ss_res, 'tss': ss_tot}

    if X.ndim == 2 and X.shape[1] == 1:
        metrics['corr'] = np.corrcoef(X.ravel(), y_true)[0,1]
    return metrics

def classification_metrics(y_true, y_pred, y_prob=None):
    metrics = {'accuracy': accuracy_score(y_true, y_pred),
               'precision': precision_score(y_true, y_pred, zero_division=0),
               'recall': recall_score(y_true, y_pred, zero_division=0),
               'f1': f1_score(y_true, y_pred, zero_division=0)}
    if y_prob is not None:
        metrics['roc_auc'] = roc_auc_score(y_true, y_prob)
    return metrics

###################### MODELING FUNCTIONS ####################################
def linear_regression(X, y):
    aligned = align_xy(X, y)
    if aligned is None:
        return None
    X, y = aligned
    if len(y) < 3:
        return None
    model = LinearRegression()
    model.fit(X, y)
    y_pred = model.predict(X)
    metrics = regression_metrics(X, y, y_pred, n_features=X.shape[1])
    return {'model': model, 'y_pred': y_pred, 'metrics': metrics}

def logistic_regression(X, y, threshold=None, regression=False):
    aligned = align_xy(X, y)
    if aligned is None:
        return None
    X, y = aligned
    y_numeric = y.copy()
    best_threshold = threshold
    regression_fit = None
    if not np.issubdtype(y.dtype, np.bool_):
        if threshold is None:
            thresholds = np.quantile(y, np.linspace(0.1, 0.9, 25))
            thresholds = np.unique(thresholds)
            best_score = -np.inf
            for thr in thresholds:
                y_class = y >= thr
                if np.unique(y_class).size < 2:
                    continue
                model = LogisticRegression()
                model.fit(X, y_class)
                y_pred = model.predict(X)
                score = f1_score(y_class, y_pred, zero_division=0)
                if score > best_score:
                    best_score = score
                    best_threshold = thr
        if best_threshold is None:
            return None
        y = y >= best_threshold
        if np.unique(y).size < 2:
            return None
    model = LogisticRegression()
    model.fit(X, y)
    y_pred = model.predict(X)
    y_prob = model.predict_proba(X)[:, 1]
    metrics = classification_metrics(y, y_pred, y_prob)
    if regression:
        mask = y.astype(bool)
        regression_fit = linear_regression(X[mask], y_numeric[mask] )
    return {'model': model,
            'threshold': best_threshold,
            'y_pred': y_pred,
            'y_prob': y_prob,
            'metrics': metrics,
            'regression_metrics': (None if regression_fit is None
                else regression_fit['metrics']),
            'regression_model': (None if regression_fit is None
                else regression_fit['model']),
            'regression_prediction': (None if regression_fit is None
                else regression_fit['y_pred'])}

def screen(features, targets, model='regression', regression=False):
    models = {}
    thresholds = {}
    regressions = {}
    rows = []
    for feat in features.columns:
        for targ in targets.columns:
            X = features[feat].values
            y = targets[targ].values
            if model == 'regression':
                fit = linear_regression(X, y)
            elif model == 'classification':
                fit = logistic_regression(X, y, regression=regression)
                thresholds[(feat, targ)] = fit['threshold']
                regressions[(feat, targ)] = {'model': fit['regression_model'],
                                             'metrics': fit['regression_metrics'],
                                             'prediction': fit['regression_prediction']}
            if fit is None:
                continue
            rows.append({'feature': feat, 'target': targ, **fit['metrics']})
            models[(feat, targ)] = fit['model']
    if not rows:
        raise ValueError('No valid models found.')
    metrics = pd.DataFrame(rows)
    if model == 'regression':
        return RegressionResults(metrics, models=models, features=features, targets=targets)
    elif model == 'classification':
        return ClassificationResults(metrics, models=models, features=features, targets=targets,
                                     thresholds=thresholds, regressions=regressions)
