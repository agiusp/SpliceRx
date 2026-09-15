"""In-memory session store with a TTL sweeper.

Single-process only. For a multi-worker deployment this needs to move to Redis
or a shared disk cache keyed by session id.
"""
from __future__ import annotations

import shutil
import tempfile
import threading
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

from .gencode import Annotation

TTL_SECONDS = 2 * 60 * 60
SWEEP_INTERVAL = 5 * 60


@dataclass
class Session:
    id: str
    tmp_dir: Path
    created: float = field(default_factory=time.time)
    last_seen: float = field(default_factory=time.time)

    # pipeline state, filled stage by stage.
    #
    # A session can hold any/all of the 3 sjdat kinds at once (junction_counts,
    # rrs_scores, gene_matrix — see services.sjdat.SJDAT_KINDS), mirroring
    # SJSurv; `active_sjdat` picks which one every downstream step below
    # reads. `junctions`/`junctions_raw` are kept as the *current* view onto
    # `sjdat[active_sjdat]`/`sjdat_raw[active_sjdat]` so every function
    # written against a single matrix (gene-set matching, MAD, feature
    # building, PCA/UMAP, heatmap) needs no change — a Sjdat is a drop-in,
    # duck-typed superset of the plain dense `services.rds.Matrix` they
    # already expected.
    sjdat_raw: Dict[str, Any] = field(default_factory=dict)   # kind -> pristine, as-loaded Sjdat
    sjdat: Dict[str, Any] = field(default_factory=dict)       # kind -> Sjdat, narrowed to samples
    # that also have a clinical row (see api/routes.py's _sync_sjdat_intersection); == sjdat_raw
    # until a clinical table is loaded
    active_sjdat: Optional[str] = None
    junctions: Optional[Any] = None      # sjdat_raw/sjdat[active_sjdat] view — see above
    junctions_raw: Optional[Any] = None
    junction_gene_index: Optional[Any] = None  # services.junction_metadata.JunctionGeneIndex —
    #                                      the fast, no-GENCODE-needed gene lookup for a
    #                                      junction-level matrix (junction_counts / rrs_scores)
    clinical: Optional[Any] = None       # services.clinical.Clinical
    annotation: Optional[Annotation] = None
    species: str = "human"
    geneset: Optional[Any] = None        # services.geneset.GeneSet
    features: Optional[Any] = None       # services.features.FeatureMatrix

    def touch(self) -> None:
        self.last_seen = time.time()


class SessionStore:
    def __init__(self, root: Optional[Path] = None):
        self._root = root or Path(tempfile.gettempdir()) / "sjvc-sessions"
        self._root.mkdir(parents=True, exist_ok=True)
        self._sessions: Dict[str, Session] = {}
        self._lock = threading.Lock()

    def create(self) -> Session:
        sid = uuid.uuid4().hex
        tmp = self._root / sid
        tmp.mkdir(parents=True, exist_ok=True)
        s = Session(id=sid, tmp_dir=tmp)
        with self._lock:
            self._sessions[sid] = s
        return s

    def get(self, sid: str) -> Optional[Session]:
        with self._lock:
            s = self._sessions.get(sid)
        if s:
            s.touch()
        return s

    def drop(self, sid: str) -> None:
        with self._lock:
            s = self._sessions.pop(sid, None)
        if s:
            shutil.rmtree(s.tmp_dir, ignore_errors=True)

    def sweep(self) -> int:
        now = time.time()
        expired: List[str] = []
        with self._lock:
            for sid, s in list(self._sessions.items()):
                if now - s.last_seen > TTL_SECONDS:
                    expired.append(sid)
                    self._sessions.pop(sid, None)
        for sid in expired:
            shutil.rmtree(self._root / sid, ignore_errors=True)
        return len(expired)

    def start_sweeper(self) -> threading.Thread:
        def loop() -> None:
            while True:
                time.sleep(SWEEP_INTERVAL)
                try:
                    self.sweep()
                except Exception:
                    pass

        t = threading.Thread(target=loop, name="sjvc-session-sweeper", daemon=True)
        t.start()
        return t


store = SessionStore()
