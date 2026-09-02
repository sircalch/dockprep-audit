"""Command-line interface for DockPrep Audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .audit import audit_pdb
from .report import render_html


def main() -> None:
    parser = argparse.ArgumentParser(prog="dockprep-audit", description="Audit a PDB receptor before molecular docking.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    audit = subparsers.add_parser("audit", help="Audit one PDB receptor and write JSON plus HTML reports.")
    audit.add_argument("pdb", type=Path, help="Input receptor PDB file")
    audit.add_argument("--output", "-o", type=Path, default=Path("dockprep-audit-results"), help="Output directory")
    args = parser.parse_args()

    if args.command == "audit":
        if not args.pdb.is_file():
            parser.error(f"PDB file not found: {args.pdb}")
        report = audit_pdb(args.pdb)
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        (args.output / "report.html").write_text(render_html(report), encoding="utf-8")
        print(f"Status: {report['summary']['status']}")
        print(f"Reports: {args.output / 'report.json'} | {args.output / 'report.html'}")
        raise SystemExit(2 if report["summary"]["high"] else 0)


if __name__ == "__main__":
    main()
