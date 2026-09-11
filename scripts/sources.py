"""
sources.py — fetch reference data from PUBLIC, open analytical-chemistry databases.

The "connect to a database" layer, kept to the two corners of analytical chemistry that are
genuinely open AND API-clean (the full sourcing map — including the many manual/paywalled ones
— is `references/databases.md`):

  Crystallography Open Database (COD)  — free CIFs; feeds crystal_engine / crystal_pxrd / cocrystal
  PubChem (PUG-REST, NCBI)              — molecule identity + COMPUTED properties (MW, XLogP, SMILES…)

stdlib-only (urllib + json), network at call time. The URL builders and response parsers are
PURE (testable offline against a saved response); only the `*_search` / `*_fetch` / `*_lookup`
wrappers touch the network, are best-effort (return empty on failure), and CITE their source URL
in the result — the "look it up + cite, don't recall" discipline (SKILL.md). NB: PubChem's
XLogP etc. are COMPUTED, and pKa is NOT reliably present — for an experimental pKa (ΔpKa input)
use the sources named in `references/databases.md`.
"""
from __future__ import annotations
import json
import re
import time
import urllib.parse
import urllib.request

COD_BASE = "https://www.crystallography.net/cod"
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_HEADERS = {"User-Agent": "analytical-figures-skill/1.0 (research)"}


def _get(url, timeout=30, retries=2, backoff=0.6):
    """GET a URL with a small retry — PubChem/E-utilities are fronted by a cloud load balancer
    (Google Cloud) reached via a CNAME chain, whose DNS can transiently fail (`getaddrinfo
    failed`) from a constrained resolver; a retry clears it. Raises the last error if all fail."""
    last = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers=_HEADERS)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:               # DNS flake, timeout, transient 5xx
            last = e
            if attempt < retries:
                time.sleep(backoff * (attempt + 1))
    raise last


# ============================================================ Crystallography Open Database
def cod_search_url(formula=None, text=None, elements=None, exclude=None, spacegroup=None, **params):
    """Build a COD search URL (format=json). `formula` in Hill notation with spaces
    ('C9 H8 O4'); `elements`/`exclude` are lists → el1..elN / nel1..nelN; extra COD params pass
    through (amin, amax, doi, space_group_number…). Pure — no network."""
    q = {"format": "json"}
    if formula:
        q["formula"] = formula
    if text:
        q["text"] = text
    if spacegroup:
        q["spacegroup"] = spacegroup
    for i, el in enumerate(elements or [], 1):
        q[f"el{i}"] = el
    for i, el in enumerate(exclude or [], 1):
        q[f"nel{i}"] = el
    q.update({k: v for k, v in params.items() if v is not None})
    return f"{COD_BASE}/result?" + urllib.parse.urlencode(q)


def cod_parse_results(text):
    """Parse a COD JSON search response (a flat list of records keyed by `file` = the COD ID)
    → list of {id, formula, sg, a, b, c, vol, doi} dicts. Pure."""
    data = json.loads(text)
    rows = data if isinstance(data, list) else data.get("result", data.get("data", []))
    out = []
    for r in rows:
        cid = r.get("file") if isinstance(r, dict) else None
        if cid is None:
            continue
        formula = (r.get("formula") or "").strip().strip("-").strip() or None
        out.append({"id": str(cid), "formula": formula,
                    **{k: r.get(k) for k in ("sg", "a", "b", "c", "vol", "doi") if r.get(k) is not None}})
    return out


def cod_cif_url(cod_id):
    return f"{COD_BASE}/{int(cod_id)}.cif"


def cod_search(formula=None, text=None, timeout=30, **kw):
    """Search COD → {hits, n, source[, warning]}. Best-effort: [] on network/parse failure.

    ⚠ HARD RULE — a COD `formula=` search returns EVERY compound with that formula: polymorphs,
    redeterminations, isomers, ESTERS, salts, and entirely unrelated molecules (C5 H8 O4 is
    glutaric acid *or* methyl ethyl oxalate; C13 H18 N4 O6 is a caffeine–glutaric cocrystal *or* a
    pyrazole-pyridine oxide). **A formula match is NOT an identity match.** Never compute, plot, or
    cite a hit until you have read its compound name — use `cod_verify(id, expect)` /
    `cod_identity(id)`, or do it in one call with `cod_search_verified(name_contains=…, formula=…)`.
    When `formula` is given the result carries a `warning` to that effect."""
    url = cod_search_url(formula=formula, text=text, **kw)
    try:
        hits = cod_parse_results(_get(url, timeout))
    except Exception as e:
        return {"hits": [], "n": 0, "source": url, "error": repr(e)}
    out = {"hits": hits, "n": len(hits), "source": url}
    if formula:
        out["warning"] = ("formula matches are NOT identity — COD returns every compound of this "
                          "formula (esters, isomers, salts, unrelated molecules). Verify each hit's "
                          "compound name with cod_identity()/cod_verify() (or use "
                          "cod_search_verified()) before computing or citing it.")
    return out


def cod_fetch_cif(cod_id, path=None, timeout=60):
    """Download a COD CIF by ID → {id, cif, identity, source, path}. `identity` is the parsed
    compound name/title/moiety (`cif_identity`) so every fetch surfaces WHAT you actually got —
    check it against what you meant before use. Writes to `path` if given, so you can then
    `crystal_engine.load(replace(cfg, cif_path=path))`. Raises on a network failure (the fetch is
    the point — don't silently hand back an empty CIF)."""
    url = cod_cif_url(cod_id)
    cif = _get(url, timeout)
    if path:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(cif)
    return {"id": str(cod_id), "cif": cif, "identity": cif_identity(cif), "source": url, "path": path}


# ------------------------------------------------------ CIF identity (the anti-collision guard)
def _cif_field(cif_text, key):
    """Value of a CIF data item — handles the ';'-delimited text block, quoted, and bare forms.
    Pure. Returns '' for absent / '?' / '.'."""
    m = re.search(re.escape(key) + r"\s*\n;(.*?)\n;", cif_text, re.S)
    if not m:
        m = re.search(re.escape(key) + r"[ \t]+(?:'([^']*)'|\"([^\"]*)\"|(\S.*))", cif_text)
    if not m:
        return ""
    val = re.sub(r"\s+", " ", (next((g for g in m.groups() if g), "") or "")).strip().strip("'\"")
    return "" if val in ("?", ".") else val


def cif_identity(cif_text):
    """Parse the IDENTITY of a CIF (common/systematic name, publication title, formula moiety, doi,
    space group) → dict. PURE (no network). This is what turns a bare COD id into a *checkable
    molecule name* — the guard against formula-collision mistakes."""
    return {
        "name": _cif_field(cif_text, "_chemical_name_common")
                or _cif_field(cif_text, "_chemical_name_systematic"),
        "title": _cif_field(cif_text, "_publ_section_title"),
        "moiety": _cif_field(cif_text, "_chemical_formula_moiety"),
        "doi": _cif_field(cif_text, "_journal_paper_doi"),
        "spacegroup": _cif_field(cif_text, "_symmetry_space_group_name_H-M")
                      or _cif_field(cif_text, "_space_group_name_H-M_alt"),
    }


def identity_matches(identity, expect):
    """True iff `expect` (a case-insensitive substring) appears in an identity dict's
    name / title / moiety. The predicate behind `cod_verify` / `cod_search_verified`."""
    hay = " ".join(str(identity.get(k, "")) for k in ("name", "title", "moiety")).lower()
    return expect.lower() in hay


def cod_identity(cod_id, timeout=60):
    """Fetch a COD CIF and return its identity dict (+ id, source). **Run this before trusting any
    formula-search hit** — a COD `formula=` search returns every compound of that formula, not the
    one you meant."""
    got = cod_fetch_cif(cod_id, timeout=timeout)
    return {"id": str(cod_id), "source": got["source"], **got["identity"]}


def cod_verify(cod_id, expect, timeout=60):
    """→ (ok, identity): does COD entry `cod_id` actually contain `expect` (a name substring)?
    Use to confirm a formula hit is the molecule you think it is before computing / citing it."""
    ident = cod_identity(cod_id, timeout=timeout)
    return identity_matches(ident, expect), ident


def cod_search_verified(name_contains, formula=None, text=None, max_check=30, timeout=30, **kw):
    """The SAFE formula→structure path: search COD, then keep only hits whose CIF name / title /
    moiety contains `name_contains` (case-insensitive), so a formula collision can't slip through.
    Fetches each candidate's CIF to read its name (one network call per hit, capped at `max_check`).
    Returns {hits:[{id, formula, sg, doi, name, title, moiety}], n, checked, source[, warning]}."""
    base = cod_search(formula=formula, text=text, timeout=timeout, **kw)
    kept, checked = [], 0
    for h in base.get("hits", []):
        if checked >= max_check:
            break
        checked += 1
        try:
            ident = cod_identity(h["id"], timeout=timeout)
        except Exception:
            continue
        if identity_matches(ident, name_contains):
            kept.append({**h, **{k: ident[k] for k in ("name", "title", "moiety")}})
    out = {"hits": kept, "n": len(kept), "checked": checked, "source": base.get("source")}
    if base.get("n", 0) > checked:
        out["warning"] = (f"only the first {checked} of {base['n']} formula hits were "
                          f"identity-checked — raise max_check to cover all")
    return out


# ============================================================================ PubChem (PUG-REST)
# Defaults kept to unambiguously-stable property names. SMILES is deliberately OMITTED: PubChem
# has been renaming it ("CanonicalSMILES" vs "ConnectivitySMILES"/"SMILES"), and ONE invalid
# property name 400s the whole request — so request SMILES explicitly with the name current on
# your PubChem, e.g. pubchem_properties("aspirin", props=("ConnectivitySMILES",)).
_PUBCHEM_DEFAULT_PROPS = ("MolecularFormula", "MolecularWeight", "XLogP", "InChIKey")


def pubchem_property_url(ident, props=_PUBCHEM_DEFAULT_PROPS, namespace="name"):
    """Build a PUG-REST compound-property URL (JSON). namespace: 'name'|'cid'|'smiles'|'inchikey'.
    Pure — no network."""
    return (f"{PUBCHEM_BASE}/compound/{namespace}/{urllib.parse.quote(str(ident))}"
            f"/property/{','.join(props)}/JSON")


def pubchem_parse_properties(text):
    """Parse a PUG-REST property JSON → the first Properties record (a dict). Pure."""
    props = json.loads(text).get("PropertyTable", {}).get("Properties", [])
    return props[0] if props else {}


def pubchem_properties(ident, props=_PUBCHEM_DEFAULT_PROPS, namespace="name", timeout=30):
    """Look up COMPUTED properties for a compound by name/CID/SMILES → {properties, source}.
    Best-effort ({} on failure). Values are computed (XLogP, TPSA…) — for an EXPERIMENTAL pKa
    (ΔpKa input) use a source from `references/databases.md`, not this."""
    url = pubchem_property_url(ident, props=props, namespace=namespace)
    try:
        out = pubchem_parse_properties(_get(url, timeout))
    except Exception as e:
        return {"properties": {}, "source": url, "error": repr(e)}
    return {"properties": out, "source": url}
