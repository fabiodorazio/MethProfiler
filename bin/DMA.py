import pandas as pd
import logging
import random
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

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

    # read processed samplesheet: this is for matching tissue sources to experiment IDs
    samplesheet = pd.read_csv(samplesheet_path)

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
    
    ####
    # -> keep only matching samples
    # match new_gms with column names
    if samples_to_keep:
        gsm_new = pd.read_csv('/Users/fdorazio/Desktop/Projects/Metprofiler/assets/GSE59685_SampleBarcode_Re-analyzedGSMs_NewGSMs.txt', sep = '\t')

        # Find matching barcodes
        gsm_ids = samplesheet.loc[samplesheet['tissue_source'].isin(samples_to_keep)]['sampleId']
        # subset GSMs and link to barcodes
        barcodes = gsm_new.loc[gsm_new['NEW_GSMs'].isin(gsm_ids),'Sample_Barcode']
        # subset beta matrix
        match_barcodes = betas.columns.intersection(barcodes)
        sub_betas = betas.loc[:, match_barcodes]
        print(sub_betas)

    else:
        sub_betas = betas

    ####
    #assert all(samplesheet.index == betas.columns)

    return(sub_betas, samplesheet)

betas, samplesheet = load_meth('/Users/fdorazio/Desktop/Projects/Metprofiler/test/GSE59685_betas_sub.csv',
                #'/Users/fdorazio/Desktop/Projects/Metprofiler/assets/samplesheet_processed.csv',
                '/Users/fdorazio/Desktop/Projects/Metprofiler/test/samplesheet_processed_test.csv',
                test_mode='On',
                samples_to_keep=['whole blood'])
print('#1')
print(betas)
def meth_qc(betas, samplesheet):
    print(samplesheet)
    ####
    # -> Missing values
    sample_missing = betas.isna().mean(axis=0)
    bad_samples = sample_missing[sample_missing > 0.05].index
    print("Removing poor quality samples:")
    print(bad_samples.tolist())
    betas = betas.drop(columns=bad_samples)
    samplesheet = samplesheet.drop(index=bad_samples)
    ####

    ####
    # -> remove probes with high missing values
    probe_missing = betas.isna().mean(axis=1)
    betas = betas.loc[probe_missing <= 0.05]
    ####

    ####
    # -> remove low-variance probes which add little biological meaning
    probe_sd = betas.std(axis=1)
    beta = betas.loc[probe_sd > 0.02]
    ####
    return(beta)

print('#2')
betas = meth_qc(betas, samplesheet)
print(betas)

def meth_normalise(betas, method="zscore"):
    """
    Methylation normalisation.
    z-score or quantile
    """

    if method == "zscore":
        norm_betas = betas.apply(
            lambda x: (x - x.mean()) / x.std(),
            axis=1
        )

    elif method == "quantile":
        # quantile normalization
        ranked = np.sort(betas.values, axis=0)
        mean_ranks = ranked.mean(axis=1)

        ranks = betas.rank(method="min").astype(int) - 1

        norm_betas = betas.copy()

        for col in betas.columns:
            norm_betas[col] = ranks[col].map(
                lambda r: mean_ranks[r]
            )

    else:
        raise ValueError("method must be 'zscore' or 'quantile'")

    return norm_betas

print('#3')
norm_betas = meth_normalise(betas)
print(norm_betas)

def DMA(betas,
        samplesheet,
        group1,
        group2):
    """
    Differential methylation analysis (probe-wise t-test)
    """

    # ------------------------------------------------------------------
    # Get samples for each group
    # ------------------------------------------------------------------
    g1_samples = samplesheet[
        samplesheet["Condition"] == group1
    ].index.intersection(betas.columns)

    g2_samples = samplesheet[
        samplesheet["Condition"] == group2
    ].index.intersection(betas.columns)

    print(f"{group1}: {len(g1_samples)} samples")
    print(f"{group2}: {len(g2_samples)} samples")

    results = []

    # ------------------------------------------------------------------
    # Probe-wise testing
    # ------------------------------------------------------------------
    for probe in betas.index:

        g1 = betas.loc[probe, g1_samples].dropna()
        g2 = betas.loc[probe, g2_samples].dropna()

        if len(g1) < 2 or len(g2) < 2:
            continue

        stat, pval = stats.ttest_ind(
            g1,
            g2,
            equal_var=False
        )

        delta_beta = g1.mean() - g2.mean()

        results.append({
            "Probe": probe,
            "mean_group1": g1.mean(),
            "mean_group2": g2.mean(),
            "delta_beta": delta_beta,
            "t_stat": stat,
            "p_value": pval
        })

    results = pd.DataFrame(results)

    # ------------------------------------------------------------------
    # Multiple testing correction
    # ------------------------------------------------------------------
    results["FDR"] = multipletests(
        results["p_value"],
        method="fdr_bh"
    )[1]

    results = results.sort_values("FDR")

    return results

print('#4')
res = DMA(norm_betas, samplesheet, group1="Control", group2="Disease")
print(res)