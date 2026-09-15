import type { SjdatOption } from "../api";

interface Props {
  /** Which of the 3 sjdat matrices are loaded into this session (loaded on
   *  the Data tab, or restored after a refresh) — pick one to analyse. */
  options: SjdatOption[];
  activeSjdat: SjdatOption["kind"] | null;
  onPick: (kind: SjdatOption["kind"]) => void;
  hasJunctionMetadata: boolean;
  /** Descriptions of the rest of what's loaded — shown as chips. */
  loaded?: { clinical?: string; gencode?: string };
}

/**
 * Read-only summary of what the Data tab has loaded into this session, plus
 * the matrix picker — this tab can analyse any of the 3 sjdat matrices
 * (junction counts, RRS scores, or the gene-level matrix), same as SJSurv.
 * Everything here is chosen on the Data tab; there is no file upload.
 */
export default function UploadPanel({ options, activeSjdat, onPick, hasJunctionMetadata, loaded }: Props) {
  return (
    <>
      <div className="radio-row">
        {options.map((o) => (
          <label key={o.kind} className={o.loaded ? "" : "disabled"}>
            <input
              type="radio"
              disabled={!o.loaded}
              checked={activeSjdat === o.kind && o.loaded}
              onChange={() => onPick(o.kind)}
            />
            <span>
              <b>{o.label}</b>
              {o.loaded ? (
                <span className="mono">
                  {" "}
                  — {o.n_features.toLocaleString()} {o.kind === "gene_matrix" ? "genes" : "features"} ×{" "}
                  {o.n_samples} samples{o.sparse ? " (sparse)" : ""}
                </span>
              ) : (
                <span className="muted"> — not loaded on the Data tab</span>
              )}
              <span className="desc">{o.description}</span>
            </span>
          </label>
        ))}
      </div>

      <div className="row" style={{ marginTop: 12 }}>
        <label style={{ display: "block" }}>Junction gene annotation</label>
        {hasJunctionMetadata ? (
          <span className="chip">✓ loaded</span>
        ) : (
          <span className="muted">
            Optional — load the junction metadata file on the Data tab for fast typed-gene /
            pathway lookups against a junction-level matrix (junction counts or RRS scores)
            without needing a GENCODE reference. Not needed for the gene matrix, which is
            already gene-level.
          </span>
        )}
      </div>

      <div className="row" style={{ marginTop: 12 }}>
        <label style={{ display: "block" }}>Clinical table</label>
        {loaded?.clinical ? (
          <span className="chip">✓ {loaded.clinical}</span>
        ) : (
          <span className="muted">Loaded with the cohort on the Data tab.</span>
        )}
      </div>

      <div className="row" style={{ marginTop: 12 }}>
        <label style={{ display: "block" }}>GENCODE reference</label>
        {loaded?.gencode ? (
          <span className="chip">✓ {loaded.gencode}</span>
        ) : (
          <span className="muted">
            Optional — choose a GENCODE release on the Data tab. Needed for a junction-level
            matrix's typed-gene / pathway lookups only when no junction metadata is loaded, and
            for the gene matrix's “protein-coding genes only” MAD filter.
          </span>
        )}
      </div>
    </>
  );
}
