import { useEffect, useRef, useState } from "react";
import { api } from "../api";

interface Props {
  sessionId: string;
  disabled?: boolean;
  onSubmit: (genes: string[]) => void;
}

const MAX_GENES = 8;

/**
 * Query-gene box with autocomplete. Accepts more than one gene: each entered
 * name becomes a chip, and "Draw" plots a region spanning all of them.
 */
export default function GeneInput({ sessionId, disabled, onSubmit }: Props) {
  const [genes, setGenes] = useState<string[]>([]);
  const [value, setValue] = useState("");
  const [suggestions, setSuggestions] = useState<string[]>([]);
  const [open, setOpen] = useState(false);
  const [highlight, setHighlight] = useState(-1);
  const boxRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const q = value.trim();
    if (q.length < 1) {
      setSuggestions([]);
      return;
    }
    const t = setTimeout(async () => {
      try {
        const names = await api.suggestGenes(sessionId, q);
        setSuggestions(names.filter((n) => !genes.some((g) => g.toLowerCase() === n.toLowerCase())));
        setHighlight(-1);
        setOpen(names.length > 0);
      } catch {
        setSuggestions([]);
      }
    }, 150);
    return () => clearTimeout(t);
  }, [value, sessionId, genes]);

  useEffect(() => {
    function onDocClick(e: MouseEvent) {
      if (boxRef.current && !boxRef.current.contains(e.target as Node)) setOpen(false);
    }
    document.addEventListener("mousedown", onDocClick);
    return () => document.removeEventListener("mousedown", onDocClick);
  }, []);

  function addGene(name: string) {
    const n = name.trim();
    if (!n) return;
    setGenes((cur) =>
      cur.some((g) => g.toLowerCase() === n.toLowerCase()) || cur.length >= MAX_GENES
        ? cur
        : [...cur, n],
    );
    setValue("");
    setOpen(false);
    setSuggestions([]);
  }

  function removeGene(name: string) {
    setGenes((cur) => cur.filter((g) => g !== name));
  }

  function submit() {
    const all = value.trim()
      ? [...genes, ...(genes.some((g) => g.toLowerCase() === value.trim().toLowerCase()) ? [] : [value.trim()])]
      : genes;
    if (all.length) {
      if (value.trim()) addGene(value.trim());
      onSubmit(all);
    }
  }

  function onKeyDown(e: React.KeyboardEvent) {
    if (open && suggestions.length > 0) {
      if (e.key === "ArrowDown") {
        e.preventDefault();
        return setHighlight((h) => (h + 1) % suggestions.length);
      }
      if (e.key === "ArrowUp") {
        e.preventDefault();
        return setHighlight((h) => (h - 1 + suggestions.length) % suggestions.length);
      }
      if (e.key === "Enter" && highlight >= 0) {
        e.preventDefault();
        return addGene(suggestions[highlight]);
      }
      if (e.key === "Escape") return setOpen(false);
    }
    if (e.key === "Enter") {
      e.preventDefault();
      if (value.trim()) addGene(value.trim());
      else if (genes.length) onSubmit(genes);
    } else if (e.key === "Backspace" && !value && genes.length) {
      setGenes((cur) => cur.slice(0, -1));
    }
  }

  const canDraw = genes.length > 0 || value.trim().length > 0;

  return (
    <div style={{ display: "flex", gap: 10, alignItems: "flex-end", flexWrap: "wrap" }}>
      <div ref={boxRef} style={{ position: "relative" }}>
        <label>
          Query gene(s)
          <div
            style={{
              display: "flex",
              flexWrap: "wrap",
              gap: 4,
              alignItems: "center",
              border: "1px solid var(--border)",
              borderRadius: 6,
              background: "var(--bg)",
              padding: "3px 6px",
              minWidth: 280,
            }}
          >
            {genes.map((g) => (
              <span key={g} className="chip" style={{ background: "var(--track-blue-utr)", color: "var(--text)", borderRadius: 4, padding: "1px 4px 1px 6px" }}>
                {g}
                <button
                  type="button"
                  onClick={() => removeGene(g)}
                  style={{ background: "none", border: "none", color: "inherit", cursor: "pointer", padding: "0 2px", fontSize: 13, lineHeight: 1 }}
                  aria-label={`remove ${g}`}
                >
                  ×
                </button>
              </span>
            ))}
            <input
              type="text"
              value={value}
              placeholder={genes.length ? "add another…" : "e.g. TP53"}
              autoComplete="off"
              disabled={genes.length >= MAX_GENES}
              style={{ border: "none", outline: "none", background: "transparent", color: "var(--text)", padding: "3px 2px", flex: 1, minWidth: 120 }}
              onChange={(e) => setValue(e.target.value)}
              onFocus={() => setOpen(suggestions.length > 0)}
              onKeyDown={onKeyDown}
            />
          </div>
        </label>
        {open && suggestions.length > 0 && (
          <ul className="suggest">
            {suggestions.map((s, i) => (
              <li
                key={s}
                className={i === highlight ? "active" : ""}
                onMouseEnter={() => setHighlight(i)}
                onMouseDown={(e) => {
                  e.preventDefault();
                  addGene(s);
                }}
              >
                {s}
              </li>
            ))}
          </ul>
        )}
      </div>
      <button type="button" disabled={disabled || !canDraw} onClick={submit}>
        Draw
      </button>
    </div>
  );
}
