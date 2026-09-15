import { useEffect, useMemo, useRef, useState } from "react";

interface Props {
  samples: string[];
  value: string;
  onChange: (s: string) => void;
}

const MAX_SHOWN = 50;

/**
 * Sample picker with type-to-filter autocomplete. The full sample list is
 * already client-side (from the RDS upload), so filtering is local — start
 * typing any part of an id and click / arrow-key a match. Better than a plain
 * <select> once a cohort has hundreds of samples.
 */
export default function SampleSelect({ samples, value, onChange }: Props) {
  const [query, setQuery] = useState(value);
  const [open, setOpen] = useState(false);
  const [highlight, setHighlight] = useState(-1);
  const boxRef = useRef<HTMLDivElement>(null);

  // reflect an externally-changed selection (e.g. a fresh RDS upload)
  useEffect(() => setQuery(value), [value]);

  useEffect(() => {
    function onDocClick(e: MouseEvent) {
      if (boxRef.current && !boxRef.current.contains(e.target as Node)) {
        setOpen(false);
        setQuery(value); // drop an unfinished edit
      }
    }
    document.addEventListener("mousedown", onDocClick);
    return () => document.removeEventListener("mousedown", onDocClick);
  }, [value]);

  const matches = useMemo(() => {
    const q = query.trim().toLowerCase();
    if (!q || q === value.toLowerCase()) return samples;
    return samples.filter((s) => s.toLowerCase().includes(q));
  }, [query, samples, value]);

  const shown = matches.slice(0, MAX_SHOWN);

  function choose(s: string) {
    onChange(s);
    setQuery(s);
    setOpen(false);
    setHighlight(-1);
  }

  function onKeyDown(e: React.KeyboardEvent) {
    if (e.key === "ArrowDown") {
      e.preventDefault();
      if (!open) return setOpen(true);
      setHighlight((h) => Math.min(h + 1, shown.length - 1));
    } else if (e.key === "ArrowUp") {
      e.preventDefault();
      setHighlight((h) => Math.max(h - 1, 0));
    } else if (e.key === "Enter") {
      e.preventDefault();
      if (highlight >= 0 && highlight < shown.length) choose(shown[highlight]);
      else if (shown.length === 1) choose(shown[0]);
    } else if (e.key === "Escape") {
      setOpen(false);
      setQuery(value);
    }
  }

  return (
    <label>
      Sample
      <div ref={boxRef} style={{ position: "relative", alignSelf: "flex-start", width: "min(420px, 100%)" }}>
        <input
          type="text"
          value={query}
          placeholder="type to filter…"
          autoComplete="off"
          spellCheck={false}
          style={{ width: "100%", boxSizing: "border-box" }}
          onChange={(e) => {
            setQuery(e.target.value);
            setOpen(true);
            setHighlight(-1);
          }}
          onFocus={() => setOpen(true)}
          onKeyDown={onKeyDown}
        />
        {open && (
          <ul className="suggest">
            {shown.map((s, i) => (
              <li
                key={s}
                className={i === highlight ? "active" : ""}
                onMouseEnter={() => setHighlight(i)}
                onMouseDown={(e) => {
                  e.preventDefault();
                  choose(s);
                }}
              >
                {s}
              </li>
            ))}
            {matches.length > MAX_SHOWN && (
              <li className="muted" style={{ pointerEvents: "none" }}>
                …and {matches.length - MAX_SHOWN} more — keep typing
              </li>
            )}
            {shown.length === 0 && (
              <li className="muted" style={{ pointerEvents: "none" }}>
                no sample matches “{query.trim()}”
              </li>
            )}
          </ul>
        )}
      </div>
    </label>
  );
}
