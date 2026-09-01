import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "validate_sample_sheet.py"
SPEC = importlib.util.spec_from_file_location("validate_sample_sheet", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class SampleSheetTests(unittest.TestCase):
    def write_sheet(self, directory, rows, extra_columns=()):
        path = directory / "samples.csv"
        fields = ["sample", *extra_columns, "fastq_1", "fastq_2"]
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        return path

    def test_identity_columns_default_to_sample(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = self.write_sheet(
                directory,
                [{"sample": "Donor A", "fastq_1": "/reads/a.fastq.gz", "fastq_2": ""}],
            )
            output = directory / "validated.csv"
            MODULE.check_samplesheet(source, output)
            with output.open() as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(row["donor_id"], "Donor_A")
            self.assertEqual(row["library_id"], "Donor_A")
            self.assertEqual(row["sample"], "Donor_A_T1")

    def test_donor_assigned_to_multiple_samples_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = self.write_sheet(
                directory,
                [
                    {"sample": "sample_a", "donor_id": "donor_1", "fastq_1": "/a.fastq.gz", "fastq_2": ""},
                    {"sample": "sample_b", "donor_id": "donor_1", "fastq_1": "/b.fastq.gz", "fastq_2": ""}
                ],
                extra_columns=("donor_id",),
            )
            with self.assertRaises(AssertionError):
                MODULE.check_samplesheet(source, directory / "validated.csv")


if __name__ == "__main__":
    unittest.main()
