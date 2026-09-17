interface Props {
  /** Data already in the session, loaded on the Data tab (or restored after a
   *  page refresh) — shown as chips. */
  loaded?: { rds?: string; ann?: string; jmeta?: string };
}

/**
 * Read-only summary of what the Data tab has loaded into this session. The
 * junction matrix and sample metadata are chosen on the Data tab; there is no
 * file upload here.
 */
export default function UploadPanel({ loaded }: Props) {
  return (
    <div className="panel">
      <h2>1 · Inputs</h2>

      <div className="row">
        <label style={{ display: "block" }}>Junction matrix</label>
        {loaded?.rds ? (
          <span className="chip">✓ {loaded.rds}</span>
        ) : (
          <span className="muted">Load a cohort on the Data tab.</span>
        )}
      </div>

      <div className="row" style={{ marginTop: 12 }}>
        <label style={{ display: "block" }}>GENCODE reference</label>
        {loaded?.ann ? (
          <span className="chip">✓ {loaded.ann}</span>
        ) : (
          <span className="muted">Choose a GENCODE release on the Data tab.</span>
        )}
      </div>

      <div className="row" style={{ marginTop: 12 }}>
        <label style={{ display: "block" }}>Junction metadata (optional)</label>
        {loaded?.jmeta ? (
          <span className="chip">✓ {loaded.jmeta}</span>
        ) : (
          <span className="muted">
            Load it on the Data tab to classify arcs from the aligner's own "annotated" column.
          </span>
        )}
      </div>
    </div>
  );
}
