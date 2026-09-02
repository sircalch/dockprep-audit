"""Build a local PDB exclusion registry from the two prior public studies."""

from __future__ import annotations

import csv
from io import StringIO
from pathlib import Path
from urllib.request import urlopen

SOURCES = [
    (
        "dude_receptor_prep_audit",
        "https://raw.githubusercontent.com/sircalch/dude-receptor-prep-audit/main/data/dude_targets.csv",
        "target",
    ),
    (
        "docking_reference_audit",
        "https://raw.githubusercontent.com/sircalch/docking-reference-audit/master/data/candidates.csv",
        "case_id",
    ),
]


def main() -> None:
    output = Path("benchmark/previous-study-exclusions.csv")
    rows: list[dict[str, str]] = []
    for study, url, label_field in SOURCES:
        payload = urlopen(url, timeout=30).read().decode("utf-8-sig")
        for record in csv.DictReader(StringIO(payload)):
            pdb_id = record.get("pdb_id", "").strip().upper()
            if pdb_id:
                rows.append({"pdb_id": pdb_id, "source_study": study, "source_label": record.get(label_field, ""), "source_url": url})
    unique = {(row["pdb_id"], row["source_study"]): row for row in rows}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["pdb_id", "source_study", "source_label", "source_url"])
        writer.writeheader()
        writer.writerows(sorted(unique.values(), key=lambda row: (row["pdb_id"], row["source_study"])))
    print(f"Wrote {len(unique)} prior-study exclusions to {output}")


if __name__ == "__main__":
    main()
