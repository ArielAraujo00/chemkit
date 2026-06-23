# Data Manipulation
import pandas as pd
import numpy as np

# Algorithms
from umap import UMAP
from hdbscan import HDBSCAN

# Metrics
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist, cdist
from sklearn.manifold import trustworthiness
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score

# Hyperparameters Tuning
import optuna

def pbounds_suggest(trial, pbounds):
    params = {}
    for key, bound in pbounds.items():
        if isinstance(bound, (int, float)):
            params[key] = bound
            continue
        if len(bound) == 3:
            low, high, step = bound
        else:
            low, high = bound
            step=None
        if any(isinstance(b, float) for b in bound):
            params[key] = trial.suggest_float(key, low, high, step=step)
        else:
            params[key] = trial.suggest_int(key, low, high,
                                            step=step if step is not None else 1)
    return params

def umap_hdbscan(data, return_models=False,
                 # UMAP params
                 n_components=2, n_neighbors=15, min_dist=0.1, spread=1.0,
                 metric='euclidean', init='random', densmap=False, angular_rp_forest=False,
                 random_state=None, low_memory=True, n_jobs=1,
                 # HDBSCAN prams
                 cluster_selection_epsilon=0.0, min_cluster_size=5,
                 min_samples=None, core_dist_n_jobs=1):

    # --- UMAP ---
    reducer = UMAP(n_components=int(n_components), n_neighbors=int(n_neighbors),
                   min_dist=float(min_dist), spread=float(spread),
                   metric=metric, init=init, densmap=densmap, angular_rp_forest=angular_rp_forest,
                   random_state=int(random_state) if random_state is not None else None,
                   low_memory=low_memory, n_jobs=n_jobs)
    embedding = reducer.fit_transform(data)

    # --- HDBSCAN ---
    clusterer = HDBSCAN(cluster_selection_epsilon=float(cluster_selection_epsilon),
                        min_cluster_size=int(min_cluster_size),
                        min_samples=int(min_samples) if min_samples is not None else None,
                        core_dist_n_jobs=core_dist_n_jobs)
    labels = clusterer.fit_predict(embedding)
    if return_models:
        return embedding, labels, reducer, clusterer
    return embedding, labels


def multi_objective_eval(embedding, labels, data=None, objectives=('ss',), sample_size=2000, n_neighbors=15, klow=8, khigh=16, kpwr=2):
    scores = {}  # Store results
    
    # Creat masks
    noise = (labels == -1)
    valid = ~noise
    valid_labels = labels[valid]
    valid_embedding = embedding[valid]
    k = len(np.unique(valid_labels))

    # ----- Clustering metrics -----
    if 'ss' in objectives:
        scores['silhouette'] = -1.0 if k < 2 else silhouette_score(valid_embedding, valid_labels)
        
    if 'chs' in objectives:
        scores['calinski_harabasz'] = -1.0 if k < 2 else calinski_harabasz_score(valid_embedding, valid_labels)
        
    if 'dbs' in objectives:
        scores['davies_bouldin'] = np.inf if k < 2 else davies_bouldin_score(valid_embedding, valid_labels)

    # ----- Manifold metrics -----
    manifold_metrics = {'knn', 'local', 'global'}
    if manifold_metrics.intersection(objectives) and data is None:
        raise ValueError("'data' must be provided for manifold metrics.")
    data = np.asarray(data)
    embedding = np.asarray(embedding)
    
    # --- local structure ---
    if 'local' in objectives:
        scores['local_structure'] = trustworthiness(data, embedding, n_neighbors=n_neighbors)
    
    # --- shared subsampling ---
    if len(data) > sample_size:
        idx = np.random.choice(len(data), size=sample_size, replace=False)
        data_sub = data[idx]
        emb_sub = embedding[idx]
    
    else:
        data_sub = data
        emb_sub = embedding
    
    # --- global structure ---
    if 'global' in objectives:
        high_dist = pdist(data_sub)
        low_dist = pdist(emb_sub)
        scores['global_structure'] = spearmanr(high_dist, low_dist).correlation
    
    # --- knn preservation ---
    if 'knn' in objectives:
        nn_high = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(data)
        nn_low = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(embedding)
        high_idx = nn_high.kneighbors(return_distance=False)[:, 1:]
        low_idx = nn_low.kneighbors(return_distance=False)[:, 1:]
        overlap = np.fromiter((np.intersect1d(high, low).size / n_neighbors
                               for high, low in zip(high_idx, low_idx)),
                              dtype=float, count=len(high_idx))
        scores['knn_preservation'] = overlap.mean()
            
    # ----- Objective Oriented Metrics -----
    if 'nf' in objectives:
        scores['noise_fraction'] = noise.mean()
    if 'ns' in objectives:
        scores['noise_sum'] = noise.sum()
    
    # Clusters
    if 'k' in objectives:
        if klow <= k <= khigh:
            kb = 0.0
        elif k < klow:
            kb = (klow - k)**kpwr
        else:
            kb = (k - khigh)**kpwr
        scores['k_bounds'] = kb
    
    # Entropy
    if 'ge' in objectives:
        _, counts = np.unique(valid_labels, return_counts=True)
        probs = counts / counts.sum()
        entropy = -np.sum(probs * np.log(probs + 1e-12))
        scores['gaussian_entropy'] = entropy
    return scores

_eval_map = {'silhouette': 'maximize', 'calinski_harabasz': 'maximize',
             'davies_bouldin': 'minimize', 'local_structure': 'maximize',
             'global_structure': 'maximize', 'knn_preservation': 'maximize',
             'noise_fraction': 'minimize', 'noise_sum': 'minimize',
             'k_bounds': 'minimize', 'gaussian_entropy': 'minimize'}
_objective_map = {'ss': 'silhouette', 'chs': 'calinski_harabasz', 'dbs': 'davies_bouldin',
                  'local': 'local_structure', 'global': 'global_structure', 'knn': 'knn_preservation',
                  'nf': 'noise_fraction', 'ns': 'noise_sum', 'ge': 'gaussian_entropy', 'k': 'k_bounds'}

def optuna_optimize(data, pbounds, study_name,
                    # Main config
                    target_function=umap_hdbscan, target_params={},
                    eval_function=multi_objective_eval, eval_params={}, n_runs=1,
                    # Sampler Config
                    multivariate=True, group=True, seed=None,
                    n_startup_trials=56, constant_liar=True,
                    # Study config
                    storage='sqlite:///study.db', load_if_exists=True,
                    # Parallel config
                    n_trials=280, n_jobs=14):
    
    # ----- define objective directions -----
    selected_objectives = eval_params.get('objectives', ('ss',))
    objective_names = [_objective_map[o] for o in selected_objectives]
    direction = [_eval_map[o] for o in objective_names]
    if len(direction) == 1:
        direction = direction[0]
    
    # ----- initiate sampler and model -----
    sampler = optuna.samplers.TPESampler(seed=seed, multivariate=multivariate, group=group,
                      n_startup_trials=n_startup_trials, constant_liar=constant_liar)
    
    study = optuna.create_study(sampler=sampler, study_name=study_name,
                                directions=direction if isinstance(direction, list) else None,
                                direction=direction if isinstance(direction, str) else None,
                                storage=storage, load_if_exists=load_if_exists)

    # --- objective function ---
    def objective(trial):
        # Suggest
        params = pbounds_suggest(trial, pbounds)
        
        # Run & Eval
        if n_runs > 1:
            ensemble_list = []
            labels_list = []
            for _ in range(n_runs):
                ensemble, labels = target_function(data, **params, **target_params)
                ensemble_list.append(ensemble)
                labels_list.append(labels)
            scores = eval_function(ensemble_list, labels_list, data=data, **eval_params)
            
        else:
            ensemble, labels = target_function(data, **params, **target_params)
            scores = eval_function(ensemble, labels, data=data, **eval_params)

        # Register Evaluation
        results = [scores[name] for name in objective_names]
        return results[0] if len(results) == 1 else tuple(results)
    
    # --- run optimization (parallel) ---
    study.optimize(objective, n_trials=n_trials, n_jobs=n_jobs, show_progress_bar=True)
    return study