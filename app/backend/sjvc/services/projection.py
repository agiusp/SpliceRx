"""2D sample projection from the feature matrix — PCA or UMAP."""
from __future__ import annotations

import threading
from dataclasses import dataclass
from typing import List, Optional

import numpy as np


class ProjectionError(ValueError):
    pass


# umap-learn runs on numba, and the only numba threading layer available in
# this environment (`workqueue`) is not threadsafe: two UMAP fits running at
# once in the request threadpool crash the numba runtime and wedge the whole
# process. Serialize every UMAP fit through this lock.
_UMAP_LOCK = threading.Lock()


@dataclass
class Embedding:
    method: str
    samples: List[str]
    coords: np.ndarray              # (n_samples, k)  — k>=2
    axis_labels: List[str]          # per component, e.g. "PC1 (43%)" or "UMAP-1"
    explained_variance: Optional[List[float]] = None   # PCA only, per component


def pca(X: np.ndarray, samples: List[str], n_components: int = 10) -> Embedding:
    from sklearn.decomposition import PCA

    n = X.shape[0]
    if n < 2:
        raise ProjectionError("PCA needs at least 2 samples")
    k = int(min(n_components, n - 1, X.shape[1]))
    model = PCA(n_components=k, svd_solver="full")
    coords = model.fit_transform(X)

    # deterministic sign: make each component's largest-magnitude loading positive
    for j in range(k):
        comp = model.components_[j]
        if comp[np.argmax(np.abs(comp))] < 0:
            comp *= -1
            coords[:, j] *= -1

    evr = model.explained_variance_ratio_.tolist()
    labels = [f"PC{j + 1} ({evr[j] * 100:.0f}%)" for j in range(k)]
    return Embedding(
        method="pca", samples=samples, coords=coords,
        axis_labels=labels, explained_variance=evr,
    )


def umap(
    X: np.ndarray,
    samples: List[str],
    *,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    seed: int = 42,
) -> Embedding:
    n = X.shape[0]
    if n < 5:
        raise ProjectionError("UMAP needs at least 5 samples")
    import umap as umap_lib  # lazy: heavy import

    nn = max(2, min(int(n_neighbors), n - 1))
    reducer = umap_lib.UMAP(
        n_components=2,
        n_neighbors=nn,
        min_dist=float(min_dist),
        random_state=seed,
        metric="euclidean",
    )
    with _UMAP_LOCK:
        coords = reducer.fit_transform(X)
    return Embedding(
        method="umap", samples=samples, coords=np.asarray(coords),
        axis_labels=["UMAP-1", "UMAP-2"],
    )
