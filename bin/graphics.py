import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def plot_group_counts(samplesheet, output_plot_dir, group_column="Condition"):
    '''
    distribution of samples per group
    '''

    counts = (
        samplesheet[group_column]
        .astype(str)
        .value_counts()
    )

    plt.figure(figsize=(6, 5))

    plt.bar(
        counts.index,
        counts.values
    )

    plt.xlabel(group_column)
    plt.ylabel("Number of Samples")
    plt.title("Group Counts")

    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'{output_plot_dir}/group_distr.png')
    plt.close


def plot_sex_distribution(samplesheet, output_plot_dir, sex_column="Sex"):
    '''
    distribution of samples by Gender
    '''

    counts = (
        samplesheet[sex_column]
        .astype(str)
        .value_counts()
    )

    plt.figure(figsize=(5, 5))

    plt.bar(
        counts.index,
        counts.values
    )

    plt.xlabel("Sex")
    plt.ylabel("Count")

    plt.title("Sex Distribution")
    plt.tight_layout()
    plt.savefig(f'{output_plot_dir}/sex_distr.png')
    plt.close


def plot_age_distribution(samplesheet, output_plot_dir, age_column="age", group_column="Condition"):
    '''
    same as above for age
    '''

    plt.figure(figsize=(7, 5))

    groups = samplesheet.groupby(group_column)

    for name, df in groups:

        plt.hist(
            df[age_column].dropna(),
            bins=15,
            alpha=0.5,
            label=str(name)
        )

        plt.xlabel("Age")
        plt.ylabel("Count")

        plt.title("Age Distribution")

        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{output_plot_dir}/{name}_age_distr.png')
        plt.close


def plot_pca(betas, samplesheet, output_plot_dir, group_column="Sex"):
    '''
    Principal component analysis for betas
    '''
    samplesheet = samplesheet.copy()

    samplesheet["Sample_Barcode"] = (
        samplesheet["Sample_Barcode"]
        .astype(str)
        .str.strip()
    )

    betas.columns = (
        betas.columns
        .astype(str)
        .str.strip()
    )
    common_samples = samplesheet[
        samplesheet["Sample_Barcode"].isin(betas.columns)
    ]["Sample_Barcode"]

    # --------------------------------------------------------
    # Subset beta matrix
    # --------------------------------------------------------
    X = betas[common_samples].T

    # remove probes with missing values
    X = X.dropna(axis=1)

    # remove probes with NA
    X = X.dropna(axis=1)

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)

    metadata = (samplesheet.set_index("Sample_Barcode").loc[X.index])

    pca_df = pd.DataFrame({
        "PC1": pcs[:, 0], "PC2": pcs[:, 1], "group": metadata[group_column].values})

    plt.figure(figsize=(7, 6))
    groups = pca_df.groupby("group")

    for name, df in groups:

        plt.scatter(df["PC1"],df["PC2"],label=str(name),alpha=0.7)

    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)")
    plt.title("PCA Plot")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{output_plot_dir}/PCA.png')
    plt.close



def volcano_plot(results, output_plot_dir, fdr_threshold=0.05, delta_threshold=0.10):
    '''
    Volcano plot for differential methylated regions
    '''
    df = results.copy()

    df["minus_log10_p"] = -np.log10(df["p_value"])

    significant = (
        (df["FDR"] < fdr_threshold) &
        (np.abs(df["delta_beta"]) > delta_threshold)
    )

    plt.figure(figsize=(8, 6))
    # non-significant
    plt.scatter(
        df.loc[~significant, "delta_beta"],
        df.loc[~significant, "minus_log10_p"],
        s=8,
        alpha=0.5,
        label="Not Significant")

    # significant
    plt.scatter(
        df.loc[significant, "delta_beta"],
        df.loc[significant, "minus_log10_p"],
        s=8,
        alpha=0.8,
        label="Significant")
    # thrsholds
    plt.axvline(delta_threshold,linestyle="--")
    plt.axvline(-delta_threshold,linestyle="--")
    plt.axhline(-np.log10(0.05),linestyle="--")
    plt.xlabel("Delta Beta")
    plt.ylabel("-log10(p-value)")
    plt.title("Volcano Plot")

    plt.savefig(f'{output_plot_dir}/volcano.png')
    plt.close


def plot_top_cpg(betas, output_plot_dir, samplesheet, results, probe=None, group_column="Condition"):

    if probe is None:
        probe = results.iloc[0]["Probe"]

    samplesheet = samplesheet.copy()

    samplesheet["Sample_Barcode"] = (
        samplesheet["Sample_Barcode"]
        .astype(str)
        .str.strip()
    )

    betas.columns = betas.columns.astype(str).str.strip()

    common = [s for s in betas.columns if s in samplesheet["Sample_Barcode"].values]
    # Build aligned dataframe
    meta = samplesheet.set_index("Sample_Barcode").loc[common]

    tmp = pd.DataFrame({
        "beta": betas.loc[probe, common].values,
        "group": meta[group_column].values
    })

    tmp = tmp.dropna()

    plt.figure(figsize=(6, 5))

    groups = tmp.groupby("group")

    for name, vals in groups:
        plt.scatter(
            [name] * len(vals),
            vals["beta"],
            alpha=0.7
        )

    plt.ylabel("Beta Value")
    plt.title(probe)

    plt.savefig(f'{output_plot_dir}/Top_cpg.png')
    plt.close

