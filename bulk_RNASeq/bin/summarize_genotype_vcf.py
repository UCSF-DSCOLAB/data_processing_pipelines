#!/usr/bin/env python3
"""Summarize donor genotype quality and pairwise discrimination in a cohort VCF."""

import argparse
import itertools
import subprocess
from pathlib import Path


MISSING_GENOTYPES = {".", "./.", ".|."}


def run_bcftools(*args):
    result = subprocess.run(
        ["bcftools", *args],
        check=True,
        text=True,
        capture_output=True,
    )
    return result.stdout


def canonical_gt(gt):
    if gt in MISSING_GENOTYPES:
        return None
    alleles = gt.replace("|", "/").split("/")
    if any(allele == "." for allele in alleles):
        return None
    return tuple(sorted(int(allele) for allele in alleles))


def numeric(value):
    if value in {"", "."}:
        return None
    return float(value)


def summarize_sample(vcf, sample):
    text = run_bcftools(
        "query", "-s", sample, "-f", "[%GT\\t%DP\\t%GQ\\n]", str(vcf)
    )
    counts = {"total": 0, "called": 0, "hom_ref": 0, "het": 0, "hom_alt": 0}
    depths = []
    qualities = []
    for line in text.splitlines():
        if not line:
            continue
        gt_text, dp_text, gq_text = line.split("\\t")
        counts["total"] += 1
        gt = canonical_gt(gt_text)
        if gt is None:
            continue
        counts["called"] += 1
        if len(set(gt)) > 1:
            counts["het"] += 1
        elif gt[0] == 0:
            counts["hom_ref"] += 1
        else:
            counts["hom_alt"] += 1
        dp = numeric(dp_text)
        gq = numeric(gq_text)
        if dp is not None:
            depths.append(dp)
        if gq is not None:
            qualities.append(gq)
    counts["call_rate"] = counts["called"] / counts["total"] if counts["total"] else 0
    counts["mean_dp"] = sum(depths) / len(depths) if depths else 0
    counts["mean_gq"] = sum(qualities) / len(qualities) if qualities else 0
    return counts


def pairwise_counts(vcf, samples):
    fmt = "%CHROM\\t%POS[\\t%GT]\\n"
    text = run_bcftools("query", "-f", fmt, str(vcf))
    pairs = {(a, b): [0, 0] for a, b in itertools.combinations(samples, 2)}
    for line in text.splitlines():
        fields = line.split("\\t")
        genotypes = [canonical_gt(value) for value in fields[2:]]
        for left, right in itertools.combinations(range(len(samples)), 2):
            gt_left, gt_right = genotypes[left], genotypes[right]
            if gt_left is None or gt_right is None:
                continue
            pair = (samples[left], samples[right])
            pairs[pair][0] += 1
            if gt_left != gt_right:
                pairs[pair][1] += 1
    return pairs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", type=Path)
    parser.add_argument("--sample-output", type=Path, required=True)
    parser.add_argument("--pairwise-output", type=Path, required=True)
    parser.add_argument("--minimum-discordant", type=int, default=50)
    args = parser.parse_args()

    samples = [line for line in run_bcftools("query", "-l", str(args.vcf)).splitlines() if line]
    if not samples:
        raise SystemExit("No donor samples found in VCF")

    with args.sample_output.open("w") as handle:
        handle.write(
            "sample\\ttotal_sites\\tcalled_sites\\tcall_rate\\thom_ref\\thet\\t"
            "hom_alt\\tmean_dp\\tmean_gq\\n"
        )
        for sample in samples:
            summary = summarize_sample(args.vcf, sample)
            handle.write(
                f"{sample}\\t{summary['total']}\\t{summary['called']}\\t"
                f"{summary['call_rate']:.4f}\\t{summary['hom_ref']}\\t"
                f"{summary['het']}\\t{summary['hom_alt']}\\t"
                f"{summary['mean_dp']:.2f}\\t{summary['mean_gq']:.2f}\\n"
            )

    with args.pairwise_output.open("w") as handle:
        handle.write("donor_1\\tdonor_2\\tjoint_called\\tdiscordant\\tstatus\\n")
        for (left, right), (joint, discordant) in pairwise_counts(args.vcf, samples).items():
            status = "PASS" if discordant >= args.minimum_discordant else "WARN"
            handle.write(f"{left}\\t{right}\\t{joint}\\t{discordant}\\t{status}\\n")


if __name__ == "__main__":
    main()
