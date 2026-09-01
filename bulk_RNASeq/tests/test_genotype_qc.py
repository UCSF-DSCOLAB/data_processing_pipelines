import importlib.util
import unittest
from pathlib import Path
from unittest.mock import patch


SCRIPT = Path(__file__).parents[1] / "bin" / "summarize_genotype_vcf.py"
SPEC = importlib.util.spec_from_file_location("summarize_genotype_vcf", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class GenotypeQCTests(unittest.TestCase):
    def test_canonical_gt(self):
        self.assertEqual(MODULE.canonical_gt("1|0"), (0, 1))
        self.assertEqual(MODULE.canonical_gt("1/1"), (1, 1))
        self.assertIsNone(MODULE.canonical_gt("./."))
        self.assertIsNone(MODULE.canonical_gt("0/."))

    def test_numeric(self):
        self.assertEqual(MODULE.numeric("12"), 12)
        self.assertIsNone(MODULE.numeric("."))

    def test_summarize_sample_parses_bcftools_tabs(self):
        query = "0/0\t10\t30\n0/1\t12\t40\n./.\t.\t.\n"
        with patch.object(MODULE, "run_bcftools", return_value=query):
            summary = MODULE.summarize_sample(Path("cohort.vcf.gz"), "donor_1")
        self.assertEqual(summary["total"], 3)
        self.assertEqual(summary["called"], 2)
        self.assertEqual(summary["hom_ref"], 1)
        self.assertEqual(summary["het"], 1)
        self.assertEqual(summary["call_rate"], 2 / 3)
        self.assertEqual(summary["mean_dp"], 11)

    def test_pairwise_counts_parses_bcftools_tabs(self):
        query = "1\t100\t0/0\t0/1\n1\t200\t1/1\t1/1\n1\t300\t./.\t0/0\n"
        with patch.object(MODULE, "run_bcftools", return_value=query):
            pairs = MODULE.pairwise_counts(
                Path("cohort.vcf.gz"), ["donor_1", "donor_2"]
            )
        self.assertEqual(pairs[("donor_1", "donor_2")], [2, 1])


if __name__ == "__main__":
    unittest.main()
