"""Portable HTML report generation with no external dependencies."""

from __future__ import annotations

from html import escape
from typing import Any


def render_html(report: dict[str, Any]) -> str:
    rows = []
    for finding in report["findings"]:
        examples = ", ".join(finding["examples"]) or "—"
        rows.append(f"<tr><td><code>{escape(finding['code'])}</code></td><td class='{finding['severity']}'>{finding['severity'].upper()}</td><td>{escape(finding['message'])}</td><td>{finding['count']}</td><td>{escape(examples)}</td></tr>")
    body = "".join(rows) or "<tr><td colspan='5'>No audit rules were triggered.</td></tr>"
    summary = report["summary"]
    return f"""<!doctype html>
<html lang='en'><head><meta charset='utf-8'><title>DockPrep Audit report</title>
<style>body{{font-family:system-ui,sans-serif;margin:2rem;color:#17202a}}header{{border-bottom:4px solid #2574a9}}table{{border-collapse:collapse;width:100%;margin-top:1.5rem}}th,td{{border:1px solid #ccd1d1;padding:.7rem;text-align:left;vertical-align:top}}th{{background:#eaf2f8}}.high{{color:#b03a2e;font-weight:700}}.medium{{color:#b9770e;font-weight:700}}.low{{color:#2471a3;font-weight:700}}code{{font-size:.9em}}</style>
</head><body><header><h1>DockPrep Audit</h1><p>Receptor-structure preparation quality-control report</p></header>
<h2>Summary</h2><p><strong>Status:</strong> {escape(summary['status'])} &nbsp; | &nbsp; High: {summary['high']} · Medium: {summary['medium']} · Low: {summary['low']}</p>
<p><strong>Source:</strong> {escape(report['source']['path'])}<br><strong>Atom records:</strong> {report['source']['atom_records']}<br><strong>Generated:</strong> {escape(report['generated_at'])}</p>
<aside><strong>Interpretation:</strong> This report identifies decisions that should be reviewed and documented. It does not automatically repair a structure or determine biological correctness.</aside>
<h2>Findings</h2><table><thead><tr><th>Rule</th><th>Severity</th><th>Why it matters</th><th>Count</th><th>Examples</th></tr></thead><tbody>{body}</tbody></table>
</body></html>"""
