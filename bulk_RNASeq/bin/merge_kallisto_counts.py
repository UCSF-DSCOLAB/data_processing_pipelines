#!/usr/bin/env python
#2023-07-31 03:14:00.000000000 +0100

"""Provide a command line tool to merge, transform and format tabular transcriptomic quantification data."""

import os
import argparse
import logging
import sys
import glob
import pandas as pd


logger = logging.getLogger()


def merge_and_transform_counts(input_count_filepaths: list):
    """
    """
    counts_df = load_counts_data(input_count_filepaths=input_count_filepaths)
    # gene_mapper_df = pd.read_csv(gene_mapper_filepath)
    # gene_mapper = dict(zip(gene_mapper_df['initial_alias'], gene_mapper_df['name']))
    # counts_df['target_id'] = (counts_df['target_id']
    #                           .apply(lambda x: gene_mapper.get(x, 'UNKNOWN')))
    gene_counts_df = (counts_df
                     .groupby(['sample_id', 'target_id'])
                     .agg(counts=('est_counts', 'sum'))
                     .reset_index()
                     )
    gene_counts_wide_df = (gene_counts_df
                   .sort_values(by=['sample_id', 'target_id'])
                   .pivot_table(index='target_id', columns='sample_id', values='counts')
                   .rename_axis(None, axis=1)
                   .astype(int)
                   .reset_index())
    return (gene_counts_wide_df
    .set_index('target_id')
    .rename_axis(None, axis=0)
    )


def load_counts_data(input_count_filepaths: list):
    return (pd
            .concat([(pd.read_csv(fp, sep='\t')
                    .assign(sample_id=fp.split('/')[-1].split('.')[0])) 
                    for fp in input_count_filepaths.split(',')]))


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge, transform, and format transcriptomic count results.",
        epilog="Example: python merge_kallisto_counts.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "input_abundance_filepaths",
        metavar="INPUT_ABUNDANCE_FILEPATHS",
        type=str,
        help="List (comma-separated) of filepaths of transcriptome quantification results for each sample.",
    )
    # parser.add_argument(
    #     "gene_mapper_filepath",
    #     metavar="INPUT_GENE_MAP",
    #     type=str,
    #     help="Filepath (CSV) for formatting the gene IDs from NCBIs standard naming to user-friendly naming scheme.",
    # )
    parser.add_argument(
        "output_merged_filepath",
        metavar="OUTPUT_MERGED_FILEPATH",
        type=str,
        help="Filepath (TSV) for the output containing merged and formatted count data.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    for filepath in args.input_abundance_filepaths.split(','):
        if not os.path.isfile(filepath):
            logger.error(f"The given input filepath {filepath} was not found!")
            sys.exit(2)
    merged_counts_df = merge_and_transform_counts(args.input_abundance_filepaths)
    merged_counts_df.to_csv(args.output_merged_filepath, sep='\t')


if __name__ == "__main__":
    sys.exit(main())
