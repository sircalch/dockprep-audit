import tempfile
import unittest
from pathlib import Path

from dockprep_audit import audit_pdb


class AuditTests(unittest.TestCase):
    def test_detects_preparation_decisions(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            pdb = Path(folder) / "example.pdb"
            pdb.write_text(
                "ATOM      1  CA AALA A   1      11.000  12.000  13.000  0.50 20.00           C  \n"
                "ATOM      2  CA BALA A   1      11.100  12.000  13.000  0.50 20.00           C  \n"
                "HETATM    3  O   HOH A 101      10.000  10.000  10.000  1.00 20.00           O  \n"
                "HETATM    4 ZN    ZN A 201      10.000  10.000  10.000  1.00 20.00          ZN  \n",
                encoding="utf-8",
            )
            report = audit_pdb(pdb)
        codes = {finding["code"] for finding in report["findings"]}
        self.assertTrue({"ALTLOC_PRESENT", "WATERS_PRESENT", "METAL_PRESENT"} <= codes)
        self.assertEqual(report["summary"]["status"], "review_required")

    def test_site_local_findings_require_ligand(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            pdb = Path(folder) / "example.pdb"
            pdb.write_text(
                "ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00 20.00           C  \n"
                "ATOM      2  CA BALA A   2      20.000  20.000  20.000  0.50 20.00           C  \n"
                "HETATM    3  O   HOH A 101      11.500  10.000  10.000  1.00 20.00           O  \n"
                "HETATM    4 ZN    ZN A 201      21.000  20.000  20.000  1.00 20.00          ZN  \n"
                "HETATM    5  C1  LIG A 300      10.000  10.000  11.500  1.00 20.00           C  \n",
                encoding="utf-8",
            )
            without_ligand = audit_pdb(pdb)
            self.assertNotIn(
                "SITE_BRIDGING_WATER_PRESENT",
                {f["code"] for f in without_ligand["findings"]},
            )

            with_ligand = audit_pdb(pdb, ligand={"component_id": "LIG", "chain": "A", "resseq": "300"})
            codes = {f["code"] for f in with_ligand["findings"]}
            # HOH A101 is 1.5 A from the ligand and 1.5 A from the ATOM A1 receptor
            # atom -- within the 3.0 A bridging-water criterion on both sides.
            self.assertIn("SITE_BRIDGING_WATER_PRESENT", codes)
            # ZN A201 and the altloc at A2 are >6 A from the ligand, so no
            # site-local finding should fire for them.
            self.assertNotIn("SITE_METAL_PRESENT", codes)
            self.assertNotIn("SITE_ALTLOC_PRESENT", codes)

    def test_site_local_findings_missing_ligand_reported(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            pdb = Path(folder) / "example.pdb"
            pdb.write_text(
                "ATOM      1  CA  ALA A   1      10.000  10.000  10.000  1.00 20.00           C  \n",
                encoding="utf-8",
            )
            report = audit_pdb(pdb, ligand={"component_id": "LIG", "chain": "A", "resseq": "300"})
            codes = {f["code"] for f in report["findings"]}
            self.assertIn("SITE_LIGAND_NOT_FOUND", codes)
