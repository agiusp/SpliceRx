import { useCallback, useEffect, useState } from "react";
import SjvApp from "./sjv/App";
import SjvcApp from "./sjvc/App";
import SjsurvApp from "./sjsurv/App";
import SjlookupApp from "./sjlookup/App";
import { api as sjvApi } from "./sjv/api";
import { api as sjvcApi } from "./sjvc/api";
import { api as sjsurvApi } from "./sjsurv/api";
import { api as sjlookupApi } from "./sjlookup/api";
import DataLoad from "./dataload/DataLoad";
import "./sjv/index.css";
import "./sjvc/index.css";
import "./sjsurv/index.css";
import "./sjlookup/index.css";

type Tab = "data" | "sjv" | "sjvc" | "sjsurv" | "sjlookup";
const TABS: Tab[] = ["data", "sjv", "sjvc", "sjsurv", "sjlookup"];

const KEY = "sj.tab";

function readStored(key: string): string | null {
  try {
    return localStorage.getItem(key);
  } catch {
    return null;
  }
}
function writeStored(key: string, value: string) {
  try {
    localStorage.setItem(key, value);
  } catch {
    /* ignore */
  }
}

/**
 * Own an app's session id here (so the Data tab can load files into the same
 * session the app tab reads from) and reuse it across page reloads — falling
 * back to a fresh session if the stored one expired or the server restarted.
 * `nonce` is bumped by `reload()` to tell the app to re-hydrate from the server.
 */
function useOwnedSession(
  storageKey: string,
  create: () => Promise<string>,
  probe: (id: string) => Promise<unknown>,
) {
  const [sid, setSid] = useState<string | null>(null);
  const [nonce, setNonce] = useState(0);

  useEffect(() => {
    const adopt = (id: string) => {
      setSid(id);
      writeStored(storageKey, id);
    };
    const fresh = () => create().then(adopt).catch(() => setSid(null));
    const stored = readStored(storageKey);
    if (stored) probe(stored).then(() => adopt(stored)).catch(fresh);
    else fresh();
  }, []); // eslint-disable-line react-hooks/exhaustive-deps

  const reload = useCallback(() => setNonce((n) => n + 1), []);
  return { sid, nonce, reload };
}

function initialTab(): Tab {
  const h = location.hash.slice(1);
  if ((TABS as string[]).includes(h)) return h as Tab;
  const s = readStored(KEY);
  if (s && (TABS as string[]).includes(s)) return s as Tab;
  return "data";
}

export default function Shell() {
  const [tab, setTab] = useState<Tab>(initialTab);

  const sjv = useOwnedSession("sj.sjvSid", sjvApi.createSession, sjvApi.sessionState);
  const sjvc = useOwnedSession("sj.sjvcSid", sjvcApi.createSession, sjvcApi.sessionState);
  const sjsurv = useOwnedSession("sj.sjsurvSid", sjsurvApi.createSession, sjsurvApi.sessionState);
  const sjlookup = useOwnedSession("sj.sjlookupSid", sjlookupApi.createSession, sjlookupApi.sessionState);

  useEffect(() => {
    writeStored(KEY, tab);
    if (location.hash.slice(1) !== tab) location.hash = tab;
  }, [tab]);

  useEffect(() => {
    const onHash = () => {
      const h = location.hash.slice(1);
      if ((TABS as string[]).includes(h)) setTab(h as Tab);
    };
    window.addEventListener("hashchange", onHash);
    return () => window.removeEventListener("hashchange", onHash);
  }, []);

  const tabButton = (t: Tab, label: string, sub: string) => (
    <button
      type="button"
      className={"sj-tab" + (tab === t ? " active" : "")}
      aria-current={tab === t}
      onClick={() => setTab(t)}
    >
      {label}
      <span className="sj-tab-sub">{sub}</span>
    </button>
  );

  return (
    <div className="sj-root">
      <header className="sj-tabs">
        {tabButton("data", "Data", "load a cohort")}
        {tabButton("sjv", "Sashimi plot", "SJV")}
        {tabButton("sjvc", "2D View", "cohort samples")}
        {tabButton("sjsurv", "Survivor groups", "SJSurv")}
        {tabButton("sjlookup", "SJ Lookup", "per-junction info")}
      </header>

      {/* Every tab stays mounted so in-progress work is preserved across a switch. */}
      <div className="sjvc-scope" hidden={tab !== "data"}>
        <DataLoad
          sjvSessionId={sjv.sid}
          sjvcSessionId={sjvc.sid}
          sjsurvSessionId={sjsurv.sid}
          sjlookupSessionId={sjlookup.sid}
          onLoaded={(target) => {
            // reload unconditionally: a partial failure can still have loaded
            // some files into the target's session server-side. One "Load
            // cohort" click can populate several tabs at once now, so we no
            // longer auto-switch here — the Data tab's own summary says which
            // tabs got data, and the tab bar above is one click away.
            if (target === "sjv") sjv.reload();
            else if (target === "sjsurv") sjsurv.reload();
            else if (target === "sjlookup") sjlookup.reload();
            else sjvc.reload();
          }}
          onGencodeApplied={() => {
            sjv.reload();
            sjvc.reload();
            sjsurv.reload();
          }}
        />
      </div>
      <div className="sjv-scope" hidden={tab !== "sjv"}>
        <SjvApp sessionId={sjv.sid} reloadNonce={sjv.nonce} />
      </div>
      <div className="sjvc-scope" hidden={tab !== "sjvc"}>
        <SjvcApp sessionId={sjvc.sid} reloadNonce={sjvc.nonce} />
      </div>
      <div className="sjsurv-scope" hidden={tab !== "sjsurv"}>
        <SjsurvApp
          sessionId={sjsurv.sid}
          reloadNonce={sjsurv.nonce}
          sjvSessionId={sjv.sid}
          sjvcSessionId={sjvc.sid}
          onGroupsPushed={(target) => {
            if (target === "sjv") sjv.reload();
            else sjvc.reload();
          }}
        />
      </div>
      <div className="sjlookup-scope" hidden={tab !== "sjlookup"}>
        <SjlookupApp sessionId={sjlookup.sid} reloadNonce={sjlookup.nonce} />
      </div>
    </div>
  );
}
