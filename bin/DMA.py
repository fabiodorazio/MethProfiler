import pandas as pd
import logging
import random
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from scipy.stats import ttest_ind

import utils
import graphics

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
    gsm_new = pd.read_csv('/Users/fdorazio/Desktop/Projects/Metprofiler/assets/GSE59685_SampleBarcode_Re-analyzedGSMs_NewGSMs.txt', sep = '\t')
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

#betas, samplesheet = load_meth('/Users/fdorazio/Desktop/Projects/Metprofiler/assets/GSE59685_betas.csv',
#                '/Users/fdorazio/Desktop/Projects/Metprofiler/assets/samplesheet_processed.csv',
#                test_mode='Off',
#                samples_to_keep=['whole blood'])

#samplesheet.to_csv('/Users/fdorazio/Desktop/Projects/Metprofiler/test/samplesheet_processed.csv')

# to delete
betas = pd.read_csv('/Users/fdorazio/Desktop/Projects/Metprofiler/test/GSE59685_betas_test_blood.csv', index_col=0)
samplesheet = pd.read_csv('/Users/fdorazio/Desktop/Projects/Metprofiler/assets/samplesheet_processed.csv')
gsm_new = pd.read_csv('/Users/fdorazio/Desktop/Projects/Metprofiler/assets/GSE59685_SampleBarcode_Re-analyzedGSMs_NewGSMs.txt', sep = '\t')
ids = gsm_new[gsm_new['Sample_Barcode'].isin(betas.columns)]['NEW_GSMs']
print(ids)
samplesheet = samplesheet[samplesheet['sampleId'].isin(ids)]
samplesheet = samplesheet.merge(
        gsm_new[['NEW_GSMs', 'Sample_Barcode']],
        left_on='sampleId',
        right_on='NEW_GSMs',
        how='left')
print('#1')
print(betas)

print(samplesheet)
OUTPUT = '/Users/fdorazio/Desktop/Projects/Metprofiler/results'
output_dir = utils.check_output_dir(OUTPUT)
output_plot_dir = utils.check_output_dir(output_dir + '/Plot_outputs')
graphics.plot_age_distribution(samplesheet,output_plot_dir, age_column="age", group_column="Condition")
graphics.plot_group_counts(samplesheet, output_plot_dir, group_column="Condition")
graphics.plot_sex_distribution(samplesheet, output_plot_dir, sex_column="Sex")
graphics.plot_pca(betas, samplesheet, output_plot_dir, group_column="Sex")

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
    removed_lowvar = (probe_sd <= 0.02).sum()
    betas = betas.loc[probe_sd > 0.02]
    print(f"Removed {removed_lowvar} low-variance probes")

    ####
    return(betas)

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


def DMA(betas, samplesheet, group1, group2):

    # ------------------------------------------------------------
    # Clean sample IDs
    # ------------------------------------------------------------
    samplesheet["Sample_Barcode"] = (
        samplesheet["Sample_Barcode"]
        .astype(str)
        .str.strip()
    )

    betas.columns = betas.columns.astype(str).str.strip()

    # ------------------------------------------------------------
    # Sample selection
    # ------------------------------------------------------------
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

    # ------------------------------------------------------------
    # Convert to matrices
    # ------------------------------------------------------------
    g1_matrix = betas[g1_samples].to_numpy(dtype=float)
    g2_matrix = betas[g2_samples].to_numpy(dtype=float)

    # ------------------------------------------------------------
    # Means
    # ------------------------------------------------------------
    mean1 = np.nanmean(g1_matrix, axis=1)
    mean2 = np.nanmean(g2_matrix, axis=1)

    delta_beta = mean1 - mean2

    # ------------------------------------------------------------
    # Vectorized t-test
    # ------------------------------------------------------------
    stat, pval = ttest_ind(
        g1_matrix,
        g2_matrix,
        axis=1,
        equal_var=False,
        nan_policy="omit"
    )

    # ------------------------------------------------------------
    # Results dataframe
    # ------------------------------------------------------------
    results = pd.DataFrame({
        "Probe": betas.index,
        "mean_group1": mean1,
        "mean_group2": mean2,
        "delta_beta": delta_beta,
        "t_stat": stat,
        "p_value": pval
    })

    # ------------------------------------------------------------
    # Remove invalid rows
    # ------------------------------------------------------------
    results = results.dropna(subset=["p_value"])

    # ------------------------------------------------------------
    # FDR correction
    # ------------------------------------------------------------
    results["FDR"] = multipletests(
        results["p_value"],
        method="fdr_bh"
    )[1]

    results = results.sort_values("FDR")


    return results

print('#4')
res = DMA(norm_betas, samplesheet, group1="Control", group2="Disease")
print(res)

def DMA_linear_model(
    betas,
    samplesheet,
    group_column="Condition",
    group1="Disease",
    group2="Control",
    covariates=["age", "Sex"],
    annotation=None
    ):
    # ------------------------------------------------------------
    # Clean sample IDs
    # ------------------------------------------------------------
    samplesheet = samplesheet.copy()

    samplesheet["Sample_Barcode"] = (
        samplesheet["Sample_Barcode"]
        .astype(str)
        .str.strip()
    )

    betas.columns = betas.columns.astype(str).str.strip()

    # ------------------------------------------------------------
    # Keep only samples present in betas
    # ------------------------------------------------------------
    samplesheet = samplesheet[
        samplesheet["Sample_Barcode"].isin(betas.columns)
    ]
    # ------------------------------------------------------------
    # Keep only groups of interest
    # ------------------------------------------------------------
    samplesheet = samplesheet[
        samplesheet[group_column].isin([group1, group2])
    ]

    print(samplesheet[group_column].value_counts())

    # ------------------------------------------------------------
    # Encode group variable
    # ------------------------------------------------------------
    samplesheet["group_binary"] = (
        samplesheet[group_column] == group1
    ).astype(int)

    # ------------------------------------------------------------
    # Encode sex
    # ------------------------------------------------------------
    if "Sex" in covariates:
        samplesheet["Sex"] = (
            samplesheet["Sex"]
            .astype(str)
            .str.upper()
        )

        samplesheet["Sex"] = samplesheet["Sex"].map({
            "MALE": 1,
            "FEMALE": 0
        })
        # ------------------------------------------------------------
    # Build design matrix
    # ------------------------------------------------------------
    design_cols = ["group_binary"] + covariates

    design = samplesheet[
        ["Sample_Barcode"] + design_cols
    ].dropna()

    # ------------------------------------------------------------
    # Match beta matrix order
    # ------------------------------------------------------------
    common_samples = design["Sample_Barcode"]

    design = design.set_index("Sample_Barcode")

    X = sm.add_constant(design)

    beta_matrix = betas[common_samples]

    print(f"Samples used: {beta_matrix.shape[1]}")
    print(f"Probes used: {beta_matrix.shape[0]}")

    # ============================================================
    # Probe-wise regression
    # ============================================================
    results = []

    for probe in beta_matrix.index:
        y = beta_matrix.loc[probe]

        # --------------------------------------------------------
        # Remove missing values
        # --------------------------------------------------------
        valid = y.notna()

        y_valid = y[valid]
        X_valid = X.loc[valid]

        if len(y_valid) < 5:
            continue

        try:
            model = sm.OLS(y_valid, X_valid)
            fit = model.fit()

            coef = fit.params["group_binary"]
            se = fit.bse["group_binary"]
            tval = fit.tvalues["group_binary"]
            pval = fit.pvalues["group_binary"]

            mean_group1 = np.nanmean(
                y_valid[
                    X_valid["group_binary"] == 1
                ]
            )

            mean_group2 = np.nanmean(
                y_valid[
                    X_valid["group_binary"] == 0
                ]
            )

            results.append({
                "Probe": probe,
                "mean_group1": mean_group1,
                "mean_group2": mean_group2,
                "delta_beta": mean_group1 - mean_group2,
                "effect_size": coef,
                "std_error": se,
                "t_stat": tval,
                "p_value": pval
            })

        except Exception:
            continue


        # ============================================================
    # Build results dataframe
    # ============================================================
    results = pd.DataFrame(results)

    # ------------------------------------------------------------
    # Multiple testing correction
    # ------------------------------------------------------------
    results["FDR"] = multipletests(
        results["p_value"],
        method="fdr_bh"
    )[1]

    # ------------------------------------------------------------
    # Rank by FDR
    # ------------------------------------------------------------
    results = results.sort_values("FDR")

    return results


graphics.plot_top_cpg(betas, output_plot_dir, samplesheet, res, probe=None, group_column="Condition")
graphics.volcano_plot(res,output_plot_dir, fdr_threshold=0.05, delta_threshold=0.10)
