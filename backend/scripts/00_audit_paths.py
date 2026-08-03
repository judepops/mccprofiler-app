#!/usr/bin/env python3
"""P0 — resolve and check every input the store builder will read.

Run this before build_store.py. Paths in the registry come from the handoffs and
have not all been verified; several artefacts moved during the July genome-wide
work. This reports what is actually on disk, searches for anything missing, and
exits non-zero if a required artefact cannot be found.

Pure stdlib. Any python3 runs it.

    python3 backend/scripts/00_audit_paths.py
    python3 backend/scripts/00_audit_paths.py --json
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app.paths import (  # noqa: E402
    ARTEFACTS,
    BLACKLIST,
    VIEWPOINT_DIR,
    search_paths,
)

# Deliberately generous. A cap of 5 silently hid the correct location of
# final_oligo_list.txt behind five sibling panel dirs on the first run — the
# audit reported "no good candidate" when the answer was there. If a name is
# ambiguous the operator needs to see all of it, not a truncated sample.
MAX_SEARCH_HITS = 40


def human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f}{unit}" if unit == "B" else f"{n:.1f}{unit}"
        n /= 1024.0
    return f"{n:.1f}GB"


def find_by_name(name: str) -> list[Path]:
    """Glob the known roots for a filename, so a moved artefact is reported not guessed."""
    hits: list[Path] = []
    for root in search_paths():
        if not root.is_dir():
            continue
        try:
            for hit in root.rglob(name):
                hits.append(hit)
                if len(hits) >= MAX_SEARCH_HITS:
                    return hits
        except (PermissionError, OSError):
            continue
    return hits


def audit() -> dict:
    rows: list[dict] = []
    missing_required: list[str] = []

    for art in ARTEFACTS:
        path = art.path
        found = path is not None and path.exists()
        suggestions = [] if found else [str(p) for p in find_by_name(path.name)]

        row = {
            "key": art.key,
            "phase": art.phase,
            "required": art.required,
            "what": art.what,
            "path": str(path) if path else None,
            "found": found,
            "size": None,
            "mtime": None,
            "suggestions": suggestions,
            "notes": art.notes,
        }

        if found and path.is_file():
            st = path.stat()
            row["size"] = st.st_size
            row["mtime"] = datetime.fromtimestamp(st.st_mtime, timezone.utc).strftime("%Y-%m-%d")

        if not found and art.required:
            missing_required.append(art.key)

        rows.append(row)

    # blacklisted files: report if present so their existence is visible, not silent
    blacklisted: list[dict] = []
    for name, why in BLACKLIST.items():
        for hit in find_by_name(name):
            blacklisted.append({"path": str(hit), "reason": why})

    viewpoints = sorted(str(p) for p in VIEWPOINT_DIR.glob("*viewpoints_final.tsv")) if VIEWPOINT_DIR.is_dir() else []

    return {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "artefacts": rows,
        "missing_required": missing_required,
        "blacklisted_present": blacklisted,
        "viewpoint_files": viewpoints,
    }


def render(result: dict) -> None:
    rows = result["artefacts"]
    width = max(len(r["key"]) for r in rows) + 2

    print("\nmccprofiler-app — P0 path audit")
    print(f"generated {result['generated']}\n")

    current_phase = None
    for r in sorted(rows, key=lambda x: (int(x["phase"][1:]), x["key"])):
        if r["phase"] != current_phase:
            current_phase = r["phase"]
            print(f"  [{current_phase}]")
        mark = "ok  " if r["found"] else ("MISS" if r["required"] else "miss")
        size = human(r["size"]) if r["size"] else "-"
        date = r["mtime"] or "-"
        print(f"    {mark}  {r['key']:<{width}} {size:>8}  {date}")
        if not r["found"]:
            print(f"          expected: {r['path']}")
            for s in r["suggestions"]:
                print(f"          found at: {s}")
            if not r["suggestions"]:
                print("          no candidates found by name search")

    vp = result["viewpoint_files"]
    print(f"\n  viewpoint BEDs: {len(vp)} found")
    for p in vp:
        print(f"    {Path(p).name}")

    if result["blacklisted_present"]:
        print("\n  BLACKLISTED FILES PRESENT (build_store will refuse to read these):")
        for b in result["blacklisted_present"]:
            print(f"    {b['path']}")
            print(f"      {b['reason']}")

    missing = result["missing_required"]
    print()
    if missing:
        print(f"  FAIL — {len(missing)} required artefact(s) missing: {', '.join(missing)}")
    else:
        optional_missing = [r["key"] for r in rows if not r["found"]]
        note = f" ({len(optional_missing)} optional missing)" if optional_missing else ""
        print(f"  PASS — all required artefacts resolved{note}")
    print()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", action="store_true", help="emit JSON instead of a table")
    args = ap.parse_args()

    result = audit()

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        render(result)

    return 1 if result["missing_required"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
