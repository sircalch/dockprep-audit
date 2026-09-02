import math
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from smoke_redock_case import (  # noqa: E402
    TIP3P_HOH_ANGLE_DEG,
    TIP3P_OH_LENGTH_A,
    _oriented_hydrogen_positions,
    append_oriented_bridging_waters_to_pdbqt,
    read_all_atom_coords_pdbqt,
)


def _angle_deg(a, vertex, b):
    u = tuple(a[i] - vertex[i] for i in range(3))
    v = tuple(b[i] - vertex[i] for i in range(3))
    dot = sum(x * y for x, y in zip(u, v))
    nu = math.dist(u, (0, 0, 0))
    nv = math.dist(v, (0, 0, 0))
    return math.degrees(math.acos(dot / (nu * nv)))


class OrientedWaterGeometryTests(unittest.TestCase):
    """F9 (PROJECT-ROADMAP.md section 19): the richer, hydrogen-bond-oriented
    TIP3P water representation used for the robustness check in Section 3.5.
    """

    def test_bond_length_and_angle_are_idealized_tip3p(self) -> None:
        oxygen = (0.0, 0.0, 0.0)
        target_a = (3.0, 0.0, 0.0)
        target_b = (0.0, 3.0, 0.0)
        h1, h2 = _oriented_hydrogen_positions(oxygen, target_a, target_b)
        self.assertAlmostEqual(math.dist(oxygen, h1), TIP3P_OH_LENGTH_A, places=6)
        self.assertAlmostEqual(math.dist(oxygen, h2), TIP3P_OH_LENGTH_A, places=6)
        self.assertAlmostEqual(_angle_deg(h1, oxygen, h2), TIP3P_HOH_ANGLE_DEG, places=6)

    def test_hydrogens_lean_toward_their_bridging_contacts(self) -> None:
        # h1 should be closer to target_a than to target_b, and vice versa for h2 --
        # otherwise the "orientation" would be meaningless (e.g. swapped or averaged
        # away rather than each H leaning toward the contact that motivated it).
        oxygen = (0.0, 0.0, 0.0)
        target_a = (3.0, 0.0, 0.0)
        target_b = (0.0, 3.0, 0.0)
        h1, h2 = _oriented_hydrogen_positions(oxygen, target_a, target_b)
        self.assertLess(math.dist(h1, target_a), math.dist(h1, target_b))
        self.assertLess(math.dist(h2, target_b), math.dist(h2, target_a))

    def test_degenerate_colinear_targets_do_not_collapse_the_hydrogens(self) -> None:
        # Both bridging contacts in the same direction from the oxygen (an edge
        # case a real bridging water can hit) must not produce h1 == h2.
        oxygen = (0.0, 0.0, 0.0)
        same_direction_target = (5.0, 0.0, 0.0)
        h1, h2 = _oriented_hydrogen_positions(oxygen, same_direction_target, same_direction_target)
        self.assertAlmostEqual(math.dist(oxygen, h1), TIP3P_OH_LENGTH_A, places=6)
        self.assertAlmostEqual(math.dist(oxygen, h2), TIP3P_OH_LENGTH_A, places=6)
        self.assertAlmostEqual(_angle_deg(h1, oxygen, h2), TIP3P_HOH_ANGLE_DEG, places=6)
        self.assertGreater(math.dist(h1, h2), 0.5)


class OrientedWaterPdbqtFormattingTests(unittest.TestCase):
    """Regression test for the PDBQT column-alignment bug found and fixed
    2026-09-01: the H1/H2 atom-name field was 3 characters wide instead of 4,
    shifting every subsequent column left by one and making Vina reject the
    x/y/z fields ("Coordinate ... is not valid"). This reads the written
    lines back with the SAME fixed-column parser the rest of the pipeline
    (and Vina) uses, so a re-introduced misalignment fails loudly here
    instead of only during an expensive redocking run.
    """

    def test_written_water_atoms_parse_back_to_the_same_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            pdbqt = Path(folder) / "receptor.pdbqt"
            pdbqt.write_text(
                "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00    -0.400 N\n"
                "TORSDOF 0\n",
                encoding="utf-8",
            )
            bridging_waters = [{
                "chain": "A", "resseq": "457",
                "oxygen_coord": (45.359, 11.570, 8.211),
                "nearest_ligand_atom": (44.542, 11.770, 8.669),
                "nearest_protein_atom": (45.083, 11.113, 7.417),
            }]
            append_oriented_bridging_waters_to_pdbqt(pdbqt, bridging_waters)

            lines = pdbqt.read_text(encoding="utf-8").splitlines()
            all_coords = read_all_atom_coords_pdbqt(pdbqt, model=None)

        water_lines = [l for l in lines if "HOH" in l]
        self.assertEqual(len(water_lines), 3)  # O, H1, H2
        # The bug this guards against (found 2026-09-01): the H1/H2 atom-name
        # field was 1 column narrower than the O line's, shifting every column
        # after it left by one. Python's lenient float() often still parses a
        # 1-column-shifted numeric field (it strips stray whitespace), so it
        # would NOT have caught this -- Vina's strict fixed-column PDBQT
        # parser did, rejecting a coordinate with a stray trailing character.
        # Checking exact line length (and thus exact column alignment) is the
        # signal that actually distinguishes the two.
        line_lengths = {len(l) for l in water_lines}
        self.assertEqual(len(line_lengths), 1,
                          f"O/H1/H2 lines have mismatched lengths (column misalignment): "
                          f"{[len(l) for l in water_lines]}")

        self.assertEqual(len(all_coords), 4)  # 1 protein atom + O + H1 + H2
        oxygen_parsed = all_coords[1]
        for expected, actual in zip((45.359, 11.570, 8.211), oxygen_parsed):
            self.assertAlmostEqual(expected, actual, places=3)

        h1_parsed, h2_parsed = all_coords[2], all_coords[3]
        for h_coord in (h1_parsed, h2_parsed):
            self.assertAlmostEqual(math.dist(oxygen_parsed, h_coord), TIP3P_OH_LENGTH_A, places=2)


if __name__ == "__main__":
    unittest.main()
