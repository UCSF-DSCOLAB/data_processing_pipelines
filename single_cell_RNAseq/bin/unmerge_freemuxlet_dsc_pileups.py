import argparse
import gzip
import os
import pandas as pd
import shutil

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--freemuxlet_dir_tsv', required=True, help='A tsv file containing `sample`, '
                        'and optionally `barcode_suffix`. The suffix need not '
                        'contain field separator (e.g. _2) since we will automatically apply `--`.')
    params = parser.parse_args()

    fmx_ids = pd.read_table(params.freemuxlet_dir_tsv, sep="\t", index_col=None, header=0, dtype=str)
    if not {'sample'}.issubset(set(fmx_ids.columns)):
        raise RuntimeError('freemuxlet_dir_tsv must necessarily have columns named `sample`')
    if 'barcode_suffix' not in fmx_ids.columns:
        if any(['--' in s for s in fmx_ids['sample']]):
            raise RuntimeError('samples inherently contain the `--` characters. Please explicitly '
                               'specify `barcode_suffix` in your tsv file.')
        fmx_ids['barcode_suffix'] = fmx_ids['sample']
    else:
        if any([s.startswith('-') for s in fmx_ids['barcode_suffix']]):
            raise RuntimeError('Values in `barcode_suffix` cannot begin with a `-`. We recommend '
                               'not using a field-separator at all.')

    if len(fmx_ids['barcode_suffix']) != len(set(fmx_ids['barcode_suffix'])):
            raise RuntimeError('Cannot have duplicate barcode suffixes. If freemuxlet_dir_tsv did '
                               ' not contain a column named `barcode_suffix`, you have duplicated '
                               '`sample` values and need to explicitly provide `barcode_suffix`.')

    assert os.path.exists('merged.lmix')
    assert os.path.exists('merged.clust1.vcf.gz')
    assert os.path.exists('merged.clust1.samples.gz')

    print('Processing merged.lmix')
    merged_lmix = pd.read_csv( 'merged.lmix', sep='\t', header=0,
                              index_col=0)
    merged_lmix['suffix'] = merged_lmix['BARCODE'].apply(lambda x: x.split('--')[-1])
    assert set(merged_lmix['suffix'].unique()) == set(fmx_ids['barcode_suffix'])
    merged_lmix['BARCODE'] == merged_lmix['BARCODE'].apply(lambda x: '--'.join(x.split('--')[:-1]))
    for s in fmx_ids.itertuples():
        print(f'Writing {s.sample}.lmix to disk')
        temp = merged_lmix.loc[merged_lmix['suffix']==s.barcode_suffix].copy()
        suffixlen = len(s.barcode_suffix) + 2  # The plus 2 is because `--` is the field separator
        temp['BARCODE'] = temp['BARCODE'].apply(lambda x: x[:-suffixlen])
        temp.drop('suffix', axis=1, inplace=True)
        temp.reset_index(drop=True, inplace=True)
        temp.to_csv(f'{s.sample}.lmix',
                    header=True, 
                    index=True, 
                    index_label='INT_ID',
                    sep='\t')

    print('Processing merged.clust1.vcf.gz')
    print('The VCF will be shared across all libraries....')
    for s in fmx_ids.itertuples():
        shutil.copy(f'merged.clust1.vcf.gz', 
                f'{s.sample}.clust1.vcf.gz')


    print('Processing merged.clust1.samples.gz')
    merged_clust1_samples_gz = pd.read_csv( 'merged.clust1.samples.gz', 
                                           sep='\t', header=0, index_col=0)
    merged_clust1_samples_gz['suffix'] = merged_clust1_samples_gz['BARCODE'].apply(lambda x: x.split('--')[-1])
    assert set(merged_clust1_samples_gz['suffix'].unique()) == set(fmx_ids['barcode_suffix'])
    merged_clust1_samples_gz['BARCODE'] == merged_clust1_samples_gz['BARCODE'].apply(lambda x: '--'.join(x.split('--')[:-1]))
    for s in fmx_ids.itertuples():
        print(f'Writing {s.sample}.clust1.samples.gz to disk')
        temp = merged_clust1_samples_gz.loc[merged_clust1_samples_gz['suffix']==s.barcode_suffix].copy()
        suffixlen = len(s.barcode_suffix) + 2  # The plus 2 is because `--` is the field separator
        temp['BARCODE'] = temp['BARCODE'].apply(lambda x: x[:-suffixlen])
        temp.drop('suffix', axis=1, inplace=True)
        temp.reset_index(drop=True, inplace=True)
        temp.to_csv(f'{s.sample}.clust1.samples.gz',
                    header=True, 
                    index=True, 
                    index_label='INT_ID',
                    sep='\t')


if __name__ == '__main__':
    main()
