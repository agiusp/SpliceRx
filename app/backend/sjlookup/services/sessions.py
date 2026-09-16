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
from typing import Any, Dict, Optional, List

TTL_SECONDS = 2 * 60 * 60
SWEEP_INTERVAL = 5 * 60


@dataclass
class Session:
    id: str
    tmp_dir: Path
    created: float = field(default_factory=time.time)
    last_seen: float = field(default_factory=time.time)

    index: Optional[Any] = None   # services.lookup.JunctionLookupIndex

    # optional sjdat matrices (junction_counts / rrs_scores only — the two
    # junction-level ones, keyed by "chr:start-end:strand" like `index`) —
    # kind -> sjvc.services.sjdat.Sjdat. Loading these is optional; a lookup
    # works from `index` alone, they only add the per-sample table for a
    # single-junction query.
    sjdat: Dict[str, Any] = field(default_factory=dict)

    def touch(self) -> None:
        self.last_seen = time.time()


class SessionStore:
    def __init__(self, root: Optional[Path] = None):
        self._root = root or Path(tempfile.gettempdir()) / "sjlookup-sessions"
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

        t = threading.Thread(target=loop, name="sjlookup-session-sweeper", daemon=True)
        t.start()
        return t


store = SessionStore()
