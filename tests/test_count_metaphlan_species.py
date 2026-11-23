import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "count_metaphlan_species.py"
INDEX = "mpa_test_index"


class SpeciesCountsTest(unittest.TestCase):
    def run_cli(self, table):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "table.tsv"
            target = Path(tmp) / "counts.tsv"
            source.write_text(table, encoding="utf-8")
            process = subprocess.run(
                [sys.executable, str(SCRIPT), "--input", str(source),
                 "--output", str(target), "--index", INDEX],
                capture_output=True, text=True,
            )
            return process, target.read_text() if target.exists() else None

    def test_versioned_two_samples(self):
        process, result = self.run_cli(
            f"#{INDEX}\nclade_name\tcase_a_{INDEX}\tcase_b_{INDEX}\n"
            "k__Bacteria|g__Genus\t2\t5\n"
            "k__Bacteria|g__Genus|s__One\t2\t0\n"
            "k__Bacteria|g__Genus|s__One|t__SGB1\t99\t99\n"
            "k__Bacteria|g__Genus|s__Two\t3\t4\n"
        )
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(result, "# samples\ttotal features\ncase_a\t2\ncase_b\t1\n")

    def test_legacy_header_and_no_species(self):
        process, result = self.run_cli("#clade_name\tA_B\nUNCLASSIFIED\t100\n")
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(result, "# samples\ttotal features\nA_B\t0\n")

    def test_rejects_bad_tables(self):
        for table in (
            "#mpa_v202403\n",
            "feature\tsample\n",
            f"clade_name\tA\tA_{INDEX}\n",
            "clade_name\tA\tB\ns__One\t1\n",
            "clade_name\tA\ns__One\tnan\n",
            "clade_name\tA\ns__One\t-1\n",
            "clade_name\tA\ns__One\n",
        ):
            with self.subTest(table=table):
                process, result = self.run_cli(table)
                self.assertNotEqual(process.returncode, 0)
                self.assertIn("error:", process.stderr)
                self.assertIsNone(result)


if __name__ == "__main__":
    unittest.main()
