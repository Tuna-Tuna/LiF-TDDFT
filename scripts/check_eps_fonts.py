#!/usr/bin/env python3
"""Fail on undefined EPS glyphs; optionally check the EPS -> PDF -> PNG path."""

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from lif_tddft.eps_export import render_eps_for_review, validate_eps_fonts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", type=Path, nargs="+", help="EPS files or directories of EPS files")
    parser.add_argument("--render-dir", type=Path, help="Also convert EPS to PDF and render the PDF here")
    parser.add_argument("--ghostscript", help="Ghostscript executable (otherwise search PATH)")
    args = parser.parse_args(argv)
    files = []
    for path in args.files:
        files.extend(sorted(path.glob("*.eps")) if path.is_dir() else [path])
    if not files:
        parser.error("No EPS files found")
    if args.render_dir and len({p.stem for p in files}) != len(files):
        parser.error("Render files with duplicate basenames in separate commands")
    for path in files:
        result = (render_eps_for_review(path, args.render_dir, ghostscript=args.ghostscript)
                  if args.render_dir else validate_eps_fonts(path))
        print(json.dumps(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
