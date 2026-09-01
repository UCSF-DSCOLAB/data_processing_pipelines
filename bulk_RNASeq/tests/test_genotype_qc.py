import importlib.util
import unittest
from pathlib import Path


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


if __name__ == "__main__":
    unittest.main()
