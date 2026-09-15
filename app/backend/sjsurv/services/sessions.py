"""In-memory session store with a TTL sweeper (single-process only).

Mirrors ``sjvc.services.sessions`` — its own store so a cohort loaded into one
app never leaks into another.
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

from sjvc.services.sjdat import SJDAT_KINDS  # noqa: F401 (re-exported — shared with NSJCG)

TTL_SECONDS = 2 * 60 * 60
SWEEP_INTERVAL = 5 * 60


@dataclass
class Session:
    id: str
    tmp_dir: Path
    created: float = field(default_factory=time.time)
    last_seen: float = field(default_factory=time.time)

    # stage-by-stage pipeline state
    raw_metadata: Optional[Any] = None             # services.metadata.RawMetadata
    metadata: Optional[Any] = None                 # services.metadata.Metadata (stratified)
    age_bands: Optional[Any] = None                # services.metadata.AgeBands currently applied
    min_group_n: Optional[int] = None              # min_group_n currently applied
    use_histology: bool = True                     # whether Group includes Histology
    histology_map: Dict[str, str] = field(default_factory=dict)  # raw value -> merged label
    sjdat_raw: Dict[str, Any] = field(default_factory=dict)  # kind -> pristine, as-loaded Sjdat
    sjdat: Dict[str, Any] = field(default_factory=dict)   # kind -> Sjdat, narrowed to samples
    # that also have a sample-metadata row (see _sync_sjdat_intersection); == sjdat_raw
    # until metadata is loaded
    active_sjdat: Optional[str] = None             # which kind the analysis uses

    # gene-set feature selection (mirrors sjvc.services.sessions.Session — the
    # active sjdat above is a `sjvc.services.sjdat.Sjdat`, a deliberate
    # drop-in replacement for `sjvc.services.rds.Matrix`, so sjvc's geneset /
    # features / gencode / junction-metadata machinery works against it
    # unchanged; see sjvc/services/sjdat.py's docstring)
    annotation: Optional[Any] = None               # sjvc.services.gencode.Annotation
    species: str = "human"
    junction_gene_index: Optional[Any] = None      # sjvc.services.junction_metadata.JunctionGeneIndex
    geneset: Optional[Any] = None                  # sjvc.services.geneset.GeneSet
    features: Optional[Any] = None                 # sjvc.services.features.FeatureMatrix, over
    # every sample of the active matrix — narrowed to a chosen Group's
    # labelled samples only at /select-geneset time (see services/select.py)

    selection: Optional[Any] = None               # services.select.Selection
    selected_group: Optional[str] = None           # the Group the selection was made for
    labels: Optional[Dict[str, str]] = None        # sample_id -> "Good" | "Poor"
    model: Optional[Any] = None                   # services.model.TrainedModel (saved MODEL)

    def touch(self) -> None:
        self.last_seen = time.time()


class SessionStore:
    def __init__(self, root: Optional[Path] = None):
        self._root = root or Path(tempfile.gettempdir()) / "sjsurv-sessions"
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

        t = threading.Thread(target=loop, name="sjsurv-session-sweeper", daemon=True)
        t.start()
        return t


store = SessionStore()
