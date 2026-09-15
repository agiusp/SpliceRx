import { useEffect, useRef, useState } from "react";

export interface Option {
  label: string;
  sub?: string;
}

interface Props {
  placeholder?: string;
  disabled?: boolean;
  fetchOptions: (q: string) => Promise<Option[]>;
  onPick: (label: string) => void;
  minChars?: number;
}

export default function Typeahead({ placeholder, disabled, fetchOptions, onPick, minChars = 1 }: Props) {
  const [value, setValue] = useState("");
  const [opts, setOpts] = useState<Option[]>([]);
  const [open, setOpen] = useState(false);
  const [hi, setHi] = useState(-1);
  const box = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const q = value.trim();
    if (q.length < minChars) {
      setOpts([]);
      return;
    }
    const t = setTimeout(async () => {
      try {
        const o = await fetchOptions(q);
        setOpts(o);
        setHi(-1);
        setOpen(o.length > 0);
      } catch {
        setOpts([]);
      }
    }, 160);
    return () => clearTimeout(t);
  }, [value]); // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    const h = (e: MouseEvent) => {
      if (box.current && !box.current.contains(e.target as Node)) setOpen(false);
    };
    document.addEventListener("mousedown", h);
    return () => document.removeEventListener("mousedown", h);
  }, []);

  function pick(label: string) {
    onPick(label);
    setValue("");
    setOpts([]);
    setOpen(false);
  }

  return (
    <div ref={box} style={{ position: "relative", minWidth: 260 }}>
      <input
        type="text"
        style={{ width: "100%" }}
        value={value}
        placeholder={placeholder}
        disabled={disabled}
        autoComplete="off"
        onChange={(e) => setValue(e.target.value)}
        onFocus={() => setOpen(opts.length > 0)}
        onKeyDown={(e) => {
          if (!open || !opts.length) return;
          if (e.key === "ArrowDown") {
            e.preventDefault();
            setHi((h) => (h + 1) % opts.length);
          } else if (e.key === "ArrowUp") {
            e.preventDefault();
            setHi((h) => (h - 1 + opts.length) % opts.length);
          } else if (e.key === "Enter" && hi >= 0) {
            e.preventDefault();
            pick(opts[hi].label);
          } else if (e.key === "Escape") {
            setOpen(false);
          }
        }}
      />
      {open && opts.length > 0 && (
        <ul className="suggest">
          {opts.map((o, i) => (
            <li
              key={o.label}
              className={i === hi ? "active" : ""}
              onMouseEnter={() => setHi(i)}
              onMouseDown={(e) => {
                e.preventDefault();
                pick(o.label);
              }}
            >
              {o.label}
              {o.sub && <span style={{ opacity: 0.6, marginLeft: 6 }}>{o.sub}</span>}
            </li>
          ))}
        </ul>
      )}
    </div>
  );
}
