"""``Sjdat`` / ``load_sjdat`` now live in ``sjvc.services.sjdat`` — shared with
NSJCG, which also lets a user pick from the same three matrix kinds. Re-
exported here so existing imports (``from .sjdat import ...`` in this
package, and anything upstream importing ``sjsurv.services.sjdat``) keep
working unchanged."""
from __future__ import annotations

from sjvc.services.sjdat import (  # noqa: F401 (re-exported)
    SJDAT_KINDS,
    SJDAT_META,
    Sjdat,
    SjdatError,
    load_sjdat,
)
