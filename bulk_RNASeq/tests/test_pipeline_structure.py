import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]


class PipelineStructureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pipeline = (ROOT / "bulk_rna_seq.nf").read_text()
        cls.haplotype = (ROOT / "modules" / "gatk4_haplotype_caller.nf").read_text()
        cls.filters = (ROOT / "modules" / "gatk4_variant_filter.nf").read_text()
        cls.parameters = (ROOT / "config" / "parameters.config").read_text()
        cls.complete = (ROOT / "modules" / "bcftools_complete_case_filter.nf").read_text()

    def test_reference_confidence_and_joint_genotyping_are_wired(self):
        self.assertIn("--emit-ref-confidence GVCF", self.haplotype)
        self.assertIn("GATK4_COMBINEGVCFS", self.pipeline)
        self.assertIn("GATK4_GENOTYPEGVCFS", self.pipeline)
        self.assertNotIn("BCFTOOLS_MERGE_VCF(", self.pipeline)

    def test_genotype_and_site_filters_are_present(self):
        self.assertIn("--genotype-filter-expression", self.filters)
        self.assertIn("--set-filtered-genotype-to-no-call true", self.filters)
        self.assertIn("gatk_extended_site_filters", self.filters)
        self.assertIn("ReadPosRankSum", self.filters)
        self.assertIn("MQRankSum", self.filters)

    def test_rna_filter_defaults_and_complete_case_output(self):
        self.assertIn("gatk_vf_fs_filter                     = 30.0", self.parameters)
        self.assertIn("gatk_extended_site_filters            = false", self.parameters)
        self.assertIn("emit_complete_case_vcf                = true", self.parameters)
        self.assertIn("F_MISSING=0", self.complete)
        self.assertIn("BCFTOOLS_COMPLETE_CASE_FILTER", self.pipeline)


if __name__ == "__main__":
    unittest.main()
