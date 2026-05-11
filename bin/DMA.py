import pandas as pd
import logging
import random
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from scipy.stats import ttest_ind
from pathlib import Path

import utils
import graphics

logging.basicConfig(level=logging.INFO)

'''
path to beta reads is passed through parser
The csv has an header that describes the labels
> '# row1: Sample barcodes\n' '# row2: New GSMs\n' '# row3: Old (Re-analyzed) GSMs\n' '\n'

'''

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
            # get path
            test_dir = Path.cwd() / 'test'
            test_dir = utils.check_output_dir(test_dir)
            test_betas = test_dir / 'GSE59685_betas_sub.csv'

            betas.to_csv(test_betas)
    elif test_mode == 'On':
        betas = pd.read_csv(
            betas_path,
            index_col=0,
            low_memory=False
        )
    else:
        raise ValueError('Error: test_mode must be one between `On` or `Off`')
    print(betas.head())
    print(betas.shape)

    # test#1: check for missing samples. Should be flagged but not killed
    missing = set(betas.columns) - set(samplesheet.index)
    if len(missing) > 0:
        logging.warning(f"Missing metadata for samples: {missing}")
    
    ####
    # -> keep only matching samples
    # match new_gms with column names
    assets_dir = Path.cwd() / 'assets'
    asset_gsm = assets_dir / 'GSE59685_SampleBarcode_Re-analyzedGSMs_NewGSMs.txt'

    gsm_new = pd.read_csv(asset_gsm, sep = '\t')
    # merge with samplesheet
    samplesheet = samplesheet.merge(
        gsm_new[['NEW_GSMs', 'Sample_Barcode']],
        left_on='sampleId',
        right_on='NEW_GSMs',
        how='left')
    
    if samples_to_keep:

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



def meth_qc(betas, samplesheet):
    print(samplesheet)
    ####
    # -> Missing values
    sample_missing = betas.isna().mean(axis=0)
    bad_samples = sample_missing[sample_missing > 0.05].index
    print("Removing poor quality samples:")
    print(bad_samples.tolist())
    betas = betas.drop(columns=bad_samples)
    # remove from samplesheet
    samplesheet = samplesheet[~samplesheet["Sample_Barcode"].isin(bad_samples)].copy()
    ####

    ####
    # -> remove probes with high missing values
    probe_missing = betas.isna().mean(axis=1)
    removed_missing = (probe_missing > 0.05).sum()
    betas = betas.loc[probe_missing <= 0.05]
    print(f"Removed {removed_missing} probes with high missingness")
    ####

    ####
    # -> remove low-variance probes which add little biological meaning
    probe_sd = betas.std(axis=1)
    removed_lowvar = (probe_sd <= 0.005).sum()
    betas = betas.loc[probe_sd > 0.005]
    print(f"Removed {removed_lowvar} low-variance probes")
    print(betas.shape)
    ####
    return(betas)



def meth_normalise(betas, method="zscore"):
    """
    Methylation normalisation.
    z-score or quantile
    """
    print(f'Normalising using {method}')
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




def DMA(betas, samplesheet, group1, group2, output_dir='results/'):

    # Clean sample IDs
    samplesheet["Sample_Barcode"] = (
        samplesheet["Sample_Barcode"]
        .astype(str)
        .str.strip()
    )

    betas.columns = betas.columns.astype(str).str.strip()

    # Sample selection
    g1_samples = pd.Index(
        samplesheet.loc[
            samplesheet["Condition"] == group1,
            "Sample_Barcode"
        ]
    ).intersection(betas.columns)

    g2_samples = pd.Index(
        samplesheet.loc[
            samplesheet["Condition"] == group2,
            "Sample_Barcode"
        ]
    ).intersection(betas.columns)

    print(f"{group1}: {len(g1_samples)} samples")
    print(f"{group2}: {len(g2_samples)} samples")

    # Convert to matrices
    g1_matrix = betas[g1_samples].to_numpy(dtype=float)
    g2_matrix = betas[g2_samples].to_numpy(dtype=float)
    # means
    mean1 = np.nanmean(g1_matrix, axis=1)
    mean2 = np.nanmean(g2_matrix, axis=1)

    delta_beta = mean1 - mean2

    # Vectorized t-test
    stat, pval = ttest_ind(
        g1_matrix,
        g2_matrix,
        axis=1,
        equal_var=False,
        nan_policy="omit"
    )

    # output
    results = pd.DataFrame({
        "Probe": betas.index,
        "mean_group1": mean1,
        "mean_group2": mean2,
        "delta_beta": delta_beta,
        "t_stat": stat,
        "p_value": pval
    })

    # remove invalid rows
    results = results.dropna(subset=["p_value"])

    # FDR correction
    results["FDR"] = multipletests(
        results["p_value"],
        method="fdr_bh"
    )[1]

    results = results.sort_values("FDR")
    results.to_csv(f'{output_dir}/dma.csv')

    return results



def beta_to_m_values(betas, eps=1e-6):
    """
    Convert beta values to M-values for statistical testing
    """
    betas = betas.clip(eps, 1 - eps)
    return np.log2(betas / (1 - betas))


def DMA_adjusted(betas, samplesheet,
    group1,
    group2,
    output_dir="results/",
    covariates=None, sample_col="Sample_Barcode", condition_col="Condition",
    use_m_values=True,
    min_non_na_fraction=0.8,
):
    """
    Differential methylation analysis adjusted for covariates

    covariates : list or None
        Metadata columns to adjust for, e.g. ["Gender", "Age", "Batch"]

    min_non_na_fraction : float
        Keep probes with at least this fraction of non-missing samples

    """

    if covariates is None:
        covariates = []

    betas = betas.copy()
    samplesheet = samplesheet.copy()

    # clean sample IDs
    samplesheet[sample_col] = (
        samplesheet[sample_col]
        .astype(str)
        .str.strip()
    )

    betas.columns = betas.columns.astype(str).str.strip()

    # keep only group1 and group2 samples
    samplesheet = samplesheet.loc[
        samplesheet[condition_col].isin([group1, group2])
    ].copy()

    # keep samples present in beta table
    shared_samples = pd.Index(samplesheet[sample_col]).intersection(betas.columns)

    if len(shared_samples) == 0:
        raise ValueError("No matching sample IDs between betas and samplesheet.")

    samplesheet = (
        samplesheet
        .set_index(sample_col)
        .loc[shared_samples]
        .copy()
    )

    betas = betas.loc[:, shared_samples]

    # Report sample counts
    n_group1 = (samplesheet[condition_col] == group1).sum()
    n_group2 = (samplesheet[condition_col] == group2).sum()

    print(f"{group1}: {n_group1} samples")
    print(f"{group2}: {n_group2} samples")

    if n_group1 < 2 or n_group2 < 2:
        raise ValueError("Each group needs at least 2 samples.")

    # Check covariates exist
    missing_covariates = [c for c in covariates if c not in samplesheet.columns]

    if missing_covariates:
        raise ValueError(f"Missing covariates in samplesheet: {missing_covariates}")

    # Encode condition
    # group1 = 1, group2 = 0
    samplesheet["group_indicator"] = (
        samplesheet[condition_col] == group1
    ).astype(int)

    # build  matrix
    design = samplesheet[["group_indicator"] + covariates].copy()
    # One-hot encode categorical covariates
    design = pd.get_dummies(design, drop_first=True)
    # ddd intercept
    design.insert(0, "intercept", 1.0)
    # convert booleans to floats if needed
    design = design.astype(float)
    print('Design matrix')
    print(design.columns.tolist())

    # drop samples with missing covariates
    valid_samples = ~design.isna().any(axis=1)

    if valid_samples.sum() < design.shape[1] + 1:
        raise ValueError(
            "Not enough samples after removing samples with missing covariates."
        )

    design = design.loc[valid_samples]
    samplesheet = samplesheet.loc[valid_samples]
    betas = betas.loc[:, design.index]

    # Check confounding / rank deficiency
    X = design.to_numpy(dtype=float)

    rank = np.linalg.matrix_rank(X)

    if rank < X.shape[1]:
        raise ValueError(
            "Design matrix is not full rank. This usually means your group is "
            "confounded with one or more covariates."
            "For example, all cases may be male and all controls female."
        )

    # Filter probes with too much missingness
    min_non_na = int(np.ceil(betas.shape[1] * min_non_na_fraction))
    keep_probes = betas.notna().sum(axis=1) >= min_non_na

    betas = betas.loc[keep_probes].copy()

    if betas.shape[0] == 0:
        raise ValueError("No probes left after missing-value filtering.")

    # Mean beta and delta beta for interpretation
    group1_samples = samplesheet.index[samplesheet[condition_col] == group1]
    group2_samples = samplesheet.index[samplesheet[condition_col] == group2]

    mean1 = betas[group1_samples].mean(axis=1, skipna=True)
    mean2 = betas[group2_samples].mean(axis=1, skipna=True)
    delta_beta = mean1 - mean2

    # Use M-values for testing, beta values for reporting
    if use_m_values:
        test_values = beta_to_m_values(betas)
    else:
        test_values = betas.copy()


    test_values = test_values.T  
    test_values = test_values.apply(lambda x: x.fillna(x.mean()), axis=0)

    Y = test_values.to_numpy(dtype=float)  

    # Vectorized linear regression
    XtX_inv = np.linalg.inv(X.T @ X)
    beta_hat = XtX_inv @ X.T @ Y

    fitted = X @ beta_hat
    residuals = Y - fitted

    n_samples, n_predictors = X.shape
    df_resid = n_samples - n_predictors

    if df_resid <= 0:
        raise ValueError(
            "Not enough residual degrees of freedom. "
            "You have too many covariates for the number of samples."
        )

    # Standard error, t-statistic, p-value for group effect
    group_idx = list(design.columns).index("group_indicator")

    residual_variance = np.sum(residuals ** 2, axis=0) / df_resid

    se_group = np.sqrt(
        residual_variance * XtX_inv[group_idx, group_idx]
    )

    coef_group = beta_hat[group_idx, :]

    t_stat = coef_group / se_group

    pval = 2 * stats.t.sf(np.abs(t_stat), df=df_resid)


    results = pd.DataFrame({
        "Probe": betas.index,
        "mean_group1": mean1.values,
        "mean_group2": mean2.values,
        "delta_beta": delta_beta.values,
        "coef_group": coef_group,
        "t_stat": t_stat,
        "p_value": pval,
        "n_samples": n_samples,
        "df_resid": df_resid,
    })

    results = results.dropna(subset=["p_value"])

    # FDR correction
    results["FDR"] = multipletests(
        results["p_value"],
        method="fdr_bh"
    )[1]

    results = results.sort_values("FDR")
    # save top 1000 by p-value for space saving
    results.head(1000).to_csv(f'{output_dir}/dma_confounding.csv')

    return results

