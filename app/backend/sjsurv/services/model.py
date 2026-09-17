"""Classification of ``SurviverGroup`` (Good vs Poor) from the selected sjdat
features.

* ``cross_validate`` — stratified k-fold; reports out-of-fold ROC AUC plus the
  spread of per-fold AUCs and the pooled confusion matrix / accuracy. This is
  the honest estimate of predictive power.
* ``train_full`` — one model fit on every selected sample and scored on the same
  samples (resubstitution — optimistic, flagged as such), returning the
  features ranked by absolute standardised weight. The fitted model is returned
  so the session can keep it (``MODEL``).

The estimator is L2 logistic regression on ``log1p`` features, standardised —
linear so the weights are directly interpretable, ``class_weight="balanced"`` so
an uneven Good/Poor split doesn't collapse to the majority class.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from .select import Selection

# Good is the positive class.
_POS = "Good"


class ModelError(ValueError):
    pass


def _xy(sel: Selection, labels: Dict[str, str], cov_values: Optional[np.ndarray]):
    y = np.array([1 if labels[s] == _POS else 0 for s in sel.sample_ids], dtype=int)
    X = np.log1p(np.abs(sel.values.T))          # samples x molecular features
    if cov_values is not None and cov_values.size:
        X = np.hstack([X, cov_values.T])        # covariates aren't count data — no log1p
    return X, y


def _check(y: np.ndarray, n_splits: int) -> None:
    n_pos, n_neg = int(y.sum()), int((y == 0).sum())
    if n_pos < 2 or n_neg < 2:
        raise ModelError(
            f"need at least 2 Good and 2 Poor samples (have Good={n_pos}, Poor={n_neg})"
        )
    if min(n_pos, n_neg) < n_splits:
        raise ModelError(
            f"ncv={n_splits} is more than the {min(n_pos, n_neg)} samples in the smaller "
            f"class — lower ncv"
        )


def _estimator() -> LogisticRegression:
    return LogisticRegression(max_iter=2000, class_weight="balanced", C=1.0)


@dataclass
class CVResult:
    n_samples: int
    n_good: int
    n_poor: int
    n_features: int
    n_splits: int
    auc: float                      # pooled out-of-fold AUC
    fold_aucs: List[float]
    fold_auc_mean: float
    fold_auc_sd: float
    accuracy: float
    sensitivity: float              # recall for Good
    specificity: float              # recall for Poor
    confusion: List[List[int]]      # [[TN, FP], [FN, TP]]
    baseline_accuracy: float        # always-predict-majority
    messages: List[str]


def cross_validate(
    sel: Selection, labels: Dict[str, str], n_splits: int,
    *, cov_names: Sequence[str] = (), cov_values: Optional[np.ndarray] = None,
) -> CVResult:
    if n_splits < 2:
        raise ModelError("ncv (cross-validation folds) must be >= 2")
    X, y = _xy(sel, labels, cov_values)
    _check(y, n_splits)
    n_good, n_poor = int(y.sum()), int((y == 0).sum())

    cov_note = f" plus {len(cov_names)} clinical covariate(s)" if cov_names else ""
    msgs = [
        f"Building classification models to predict SurviverGroup on {len(y)} samples "
        f"(Good={n_good}, Poor={n_poor}).",
        f"{n_splits}-fold stratified cross-validation over {sel.n_features} selected "
        f"{sel.feature_noun} feature(s){cov_note}.",
    ]

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=0)
    oof = np.full(len(y), np.nan)
    fold_aucs: List[float] = []
    for k, (tr, te) in enumerate(skf.split(X, y), start=1):
        scaler = StandardScaler().fit(X[tr])
        clf = _estimator().fit(scaler.transform(X[tr]), y[tr])
        p = clf.predict_proba(scaler.transform(X[te]))[:, 1]
        oof[te] = p
        if len(np.unique(y[te])) == 2:
            fa = float(roc_auc_score(y[te], p))
            fold_aucs.append(fa)
            msgs.append(f"  fold {k}: {len(te)} test samples, AUC {fa:.3f}")
        else:
            msgs.append(f"  fold {k}: {len(te)} test samples, AUC n/a (one class)")

    auc = float(roc_auc_score(y, oof))
    pred = (oof >= 0.5).astype(int)
    cm = confusion_matrix(y, pred, labels=[0, 1])
    tn, fp, fn, tp = (int(v) for v in cm.ravel())
    acc = accuracy_score(y, pred)
    base = max(n_good, n_poor) / len(y)
    msgs.append(
        f"Pooled out-of-fold AUC {auc:.3f}; per-fold {np.mean(fold_aucs):.3f} "
        f"± {np.std(fold_aucs):.3f}." if fold_aucs else f"Pooled out-of-fold AUC {auc:.3f}."
    )
    msgs.append(
        f"Accuracy {acc:.3f} vs {base:.3f} for always predicting the larger class."
    )

    return CVResult(
        n_samples=len(y), n_good=n_good, n_poor=n_poor, n_features=sel.n_features + len(cov_names),
        n_splits=n_splits, auc=auc, fold_aucs=fold_aucs,
        fold_auc_mean=float(np.mean(fold_aucs)) if fold_aucs else float("nan"),
        fold_auc_sd=float(np.std(fold_aucs)) if fold_aucs else float("nan"),
        accuracy=float(acc),
        sensitivity=float(tp / (tp + fn)) if (tp + fn) else float("nan"),
        specificity=float(tn / (tn + fp)) if (tn + fp) else float("nan"),
        confusion=[[tn, fp], [fn, tp]],
        baseline_accuracy=float(base),
        messages=msgs,
    )


@dataclass
class FeatureWeight:
    feature: str
    weight: float                   # signed standardised coefficient
    abs_weight: float
    direction: str                  # "higher in Good" | "higher in Poor"
    mean_good: float
    mean_poor: float


@dataclass
class TrainedModel:
    sjdat_kind: str
    group: str
    n_samples: int
    n_good: int
    n_poor: int
    n_features: int
    auc_resub: float
    accuracy_resub: float
    confusion_resub: List[List[int]]
    features: List[FeatureWeight]
    # kept so the session can reuse MODEL later
    _coef: np.ndarray = None
    _intercept: float = 0.0
    _scaler_mean: np.ndarray = None
    _scaler_scale: np.ndarray = None
    _feature_ids: List[str] = None


def train_full(
    sel: Selection, labels: Dict[str, str], *, sjdat_kind: str, group: str,
    cov_names: Sequence[str] = (), cov_values: Optional[np.ndarray] = None,
) -> TrainedModel:
    X, y = _xy(sel, labels, cov_values)
    _check(y, 2)
    n_good, n_poor = int(y.sum()), int((y == 0).sum())

    scaler = StandardScaler().fit(X)
    Xs = scaler.transform(X)
    clf = _estimator().fit(Xs, y)
    p = clf.predict_proba(Xs)[:, 1]
    pred = (p >= 0.5).astype(int)
    cm = confusion_matrix(y, pred, labels=[0, 1])
    tn, fp, fn, tp = (int(v) for v in cm.ravel())

    coef = clf.coef_.ravel()
    good_mask = y == 1
    has_cov = cov_values is not None and cov_values.size
    raw = np.vstack([sel.values, cov_values]) if has_cov else sel.values   # features x samples
    feature_ids = list(sel.feature_ids) + list(cov_names) if has_cov else list(sel.feature_ids)
    mean_good = raw[:, good_mask].mean(axis=1)
    mean_poor = raw[:, ~good_mask].mean(axis=1)

    order = np.argsort(np.abs(coef))[::-1]
    feats = [
        FeatureWeight(
            feature=feature_ids[i],
            weight=float(coef[i]),
            abs_weight=float(abs(coef[i])),
            direction="higher in Good" if coef[i] >= 0 else "higher in Poor",
            mean_good=float(mean_good[i]),
            mean_poor=float(mean_poor[i]),
        )
        for i in order
    ]

    return TrainedModel(
        sjdat_kind=sjdat_kind, group=group,
        n_samples=len(y), n_good=n_good, n_poor=n_poor, n_features=len(feature_ids),
        auc_resub=float(roc_auc_score(y, p)),
        accuracy_resub=float(accuracy_score(y, pred)),
        confusion_resub=[[tn, fp], [fn, tp]],
        features=feats,
        _coef=coef, _intercept=float(clf.intercept_[0]),
        _scaler_mean=scaler.mean_, _scaler_scale=scaler.scale_,
        _feature_ids=feature_ids,
    )
