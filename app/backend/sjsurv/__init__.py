"""SJSurv — survivor-group classification from splice-junction features.

Explores whether a cohort's splice-junction features separate the samples
labelled ``Good`` / ``Poor`` in the ``SurviverGroup`` column of the sample
metadata (a ``prepTCGAdata`` / ``surv_cohort()`` output).

Pipeline (one session, filled stage by stage):

1. sample metadata with ``Group`` + ``SurviverGroup`` columns, and one or more
   splice-junction feature matrices (``sjdat``): raw junction counts, RRS
   scores, or novel-junction-counts-per-gene — all loaded on the Data tab.
2. the user picks a ``Group`` (or "all labelled samples") and which ``sjdat`` to
   use.
3. feature selection by non-zero coverage (N), magnitude (X) and MAD rank (n).
4. a cross-validated classifier (ncv folds) of ``SurviverGroup``.
5. a model trained on every selected sample, with its features ranked by weight.

Only separated from ``sjv`` / ``sjvc`` by URL prefix (``/api/sjsurv``); it reuses
their RDS readers.
"""
