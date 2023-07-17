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
    parser.add_argument('--ignore_diff_lengths_error', action='store_true')
    params = parser.parse_args()

    fmx_ids = pd.read_table(params.freemuxlet_dir_tsv, sep="\t", index_col=None, header=0)
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


    print('Processing *.var.gz')
    _x = None
    num_lines = {}
    fails = 0
    for s in fmx_ids.itertuples():
        with gzip.open(f'{s.sample}.var.gz') as iff:
            x = len(iff.readlines())
        if _x is None:
            _x = x
        fails += _x!=x
        num_lines[s.sample] = x

    if fails > 0:
        if params.ignore_diff_lengths_error:
            longest = max(num_lines, key=num_lines.get)
            longest_idx = [i for i, x in enumerate(fmx_ids['sample']==longest) if x][0]
            if len(set(list(num_lines.values()))) != 1:
                print(f'!!Found var.gz files of different lengths, but assuming this is not reference VCF-related, and using the longest one, {longest}')
            print('Writing merged.var.gz to disk')
            shutil.copy( f'{longest}.var.gz','merged.var.gz')
        else:
            print(f'Found {fails} failures. Full list below.')
            for x, y in num_lines.items():
                print(f'{x}: {y}')
            raise RuntimeError('Freemuxlet runs may have used different VCFs for processing, or it may be just a freemuxlet bug. Check if the shorter .var files are simply missing some rows from the bottom of the longer counterparts.  If so, remove the merge directory and run again with --ignore_diff_lengths_error added.')
    else:
        print('Writing merged.var.gz to disk')
        shutil.copy(f'{s.sample}.var.gz', 'merged.var.gz')

    merged_cels = None
    merged_plps = None
    print('Processing *.cel.gz and *.plp.gz concurrently')

    existing_cells = 0
    for s in fmx_ids.itertuples():
        print(f'Processing {s.sample}')
        cel = pd.read_csv( f'{s.sample}.cel.gz', sep="\t", header=0,
                          index_col=None)
        plp = pd.read_csv(f'{s.sample}.plp.gz', sep="\t", header=0,
                          index_col=None, dtype={
                              '#DROPLET_ID': int,
                              'SNP_ID': int,
                              'ALLELES': str,
                              'BASEQS': str,
                              })
        cel['#DROPLET_ID'] = cel['#DROPLET_ID'] + existing_cells
        plp['#DROPLET_ID'] = plp['#DROPLET_ID'] + existing_cells
        existing_cells += cel.shape[0]
        cel['BARCODE'] = cel['BARCODE'].apply(lambda x: f'{x}--{s.barcode_suffix}')

        if merged_cels is None:
            merged_cels = cel
            merged_plps = plp
        else:
            merged_cels = pd.concat([merged_cels, cel])
            merged_plps = pd.concat([merged_plps, plp])

    print('Writing merged.cel.gz to disk')
    merged_cels.to_csv('merged.cel.gz',
                       header=True,
                       index=False,
                       sep="\t")
    print('Writing merged.barcodes.gz to disk')
    merged_cels.to_csv('merged.barcodes.gz',
                       columns=['BARCODE'],
                       header=False,
                       index=False,
                       sep="\t")
    
    print('Sorting merged.plp.gz')
    merged_plps.sort_values(axis=0, by=['SNP_ID', '#DROPLET_ID'], inplace=True)
    print('writing merged.plp.gz to disk')
    merged_plps.to_csv('merged.plp.gz', 
                       header=True,
                       index=False,
                       sep="\t")
    print('Done.')

if __name__ == '__main__':
    main()
