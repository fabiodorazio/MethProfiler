import pandas as pd
import logging
import random

import utils

logging.basicConfig(level=logging.INFO)

'''
path to beta reads is passed through parser
The csv has an header that describes the labels
> '# row1: Sample barcodes\n' '# row2: New GSMs\n' '# row3: Old (Re-analyzed) GSMs\n' '\n'

'''
#'/Users/fdorazio/Desktop/Projects/Metprofiler/assets/GSE59685_betas.csv'
def load_meth(betas_path, samplesheet_path, sample_on='Off', test_mode='Off', samples_to_keep=None):
    '''
    read samplesheet and beta file
    extracts metadata from rows
    checks for missing samples
    returns beta table
    '''
    # read samplesheet
    import pandas as pd

    samplesheet = pd.read_csv(samplesheet_path)
    print(samplesheet)

    # set sample IDs as index
    if 'sampleId' in samplesheet.columns:
        samplesheet = samplesheet.set_index('sampleId')
    print(samplesheet)

    # read betas outside of test mode
    if test_mode == 'Off':
        meth = pd.read_csv(
            betas_path,
            skiprows=4,
            index_col=0,
            low_memory=False
        )
        # methylation matrix as numeric
        betas = meth.iloc[2:].astype(float)

        # subsample for testing
        if sample_on == 'On':
            sampled = random.sample(list(betas.columns), 20)
            betas = betas.loc[:, sampled]
            betas.to_csv(
                '/Users/fdorazio/Desktop/Projects/Metprofiler/test/GSE59685_betas_sub.csv'
            )
    elif test_mode == 'On':
        betas = pd.read_csv(
            betas_path,
            index_col=0,
            low_memory=False
        )
    else:
        raise ValueError('Error: est_mode must be one between `On` or `Off`')
    print(betas.head())
    print(betas.shape)

    # test#1: check for missing samples. Should be flagged but not killed
    missing = set(betas.columns) - set(samplesheet.index)
    if len(missing) > 0:
        logging.warning(f"Missing metadata for samples: {missing}")
    
    # keep only matching samples
    if samples_to_keep:
        samples_subset = samplesheet.loc[samples_to_keep,]
    else:
        samples_subset = samplesheet
    valid = [x for x in betas.columns if x in samples_subset['Sample_Barcode'].values]
    betas = betas.loc[:, valid]
    samplesheet = samplesheet[samplesheet['Sample_Barcode'].isin(valid)]

    assert all(samplesheet.index == betas.columns)

    return(betas)

betas = load_meth('/Users/fdorazio/Desktop/Projects/Metprofiler/test/GSE59685_betas_sub.csv',
                '/Users/fdorazio/Desktop/Projects/Metprofiler/assets/samplesheet_processed.csv',
                test_mode='On',
                samples_to_keep='GSM1068950')
print('done')
print(betas)
def meth_qc(betas):
    # Missing values
    sample_missing = betas.isna().mean(axis=0)

    bad_samples = sample_missing[sample_missing > 0.05].index

    print("Removing poor quality samples:")
    print(bad_samples.tolist())

    betas = betas.drop(columns=bad_samples)
    samplesheet = samplesheet.drop(index=bad_samples)

    # remove probes with high missing values
    probe_missing = betas.isna().mean(axis=1)
    betas = betas.loc[probe_missing <= 0.05]

    # remove low-variance probes which add little biological meaning
    probe_sd = betas.std(axis=1)

    beta = betas.loc[probe_sd > 0.02]

def meth_normalise():
    pass


def DMA(betas):
    pass
