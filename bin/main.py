import os
import argparse
import logging

import utils
import DMA
import graphics

def get_args():
    parser = argparse.ArgumentParser(prog="metprofiler")
    parser.add_argument("betas", help = "Input path to beta values")
    parser.add_argument("samplesheet", help = "Input path to samplesheet")
    parser.add_argument("-n", dest="name", default=None, help="Output file name without the format")
    parser.add_argument("-o", dest="output_dir", default="../outputs", help="Output directory path")
    parser.add_argument("--test_mode", dest="test_mode", default='Off', help="When `On`, runs the pipeline on a smaller subset of beta probes")
    parser.add_argument("--keep", dest="to_keep", default=None, help="A list of tissue sources for filtering the dataset. Default = None")
    #parser.add_argument("--normalise", dest="normalise", default=20, type=int, help="If normalisation is not needed, turns this argument off")
    parser.add_argument("--norm_method", dest="norm_method", default="zscore", type=int, help="Which normalization method is applied to betas. Either z-score or quantiles")
    #parser.add_argument("--covariates", dest="covariates", default=20, type=int, help="Number of bins for each chromosome for visualisation purposes")
    parser.add_argument("--group1", dest="group1", default="Control", type=int, help="Category of group1 comparison")
    parser.add_argument("--group2", dest="group2", default="Disease", type=int, help="Category of group2 comparison")


    args = parser.parse_args()
    beta = args.beta
    samplesheet = args.samplesheet
    name = args.name
    output = args.output_dir
    test_mode = args.test_mode
    to_keep = args.to_keep
    normalise = args.normalise
    norm_method = args.norm_method
    covariates = args.covariates
    group1 = args.group1
    group2 = args.group2

    return (beta, samplesheet, name, output, test_mode, to_keep, normalise, norm_method, covariates, group1, group2)


if __name__ == "__main__":
    BETA, SAMPLESHEET, NAME, OUTPUT, TEST_MODE, TO_KEEP, NORMALIZE, NORM_METHOD, COVARIATES, GROUP1, GROUP2 = get_args()

    if NAME is None:
        NAME = utils.get_basename(BETA)

    # checks if output exists or create one
    output_dir = utils.check_output_dir(OUTPUT)
    # creates path for beta
    output_dma = output_dir + NAME + "_DMA.csv"
    # sets location for plot output, checks if exists or creates one
    output_plot_dir = utils.check_output_dir(output_dir + '/Plot_outputs')

    #### EXECUTION

    # 1) Load and Prepare betas and samplesheet
    betas, samplesheet = DMA.load_meth(BETA,
                SAMPLESHEET,
                test_mode=TEST_MODE,
                samples_to_keep=TO_KEEP)
    
    # 2) QC
    betas = DMA.meth_qc(betas, samplesheet)

    # 3) Normalise
    if NORMALIZE:
        betas = DMA.meth_normalise(betas=betas, method=NORM_METHOD)

    # 4) Differential Analysis
    if COVARIATES:
        results = DMA.DMA_linear_model()
    else:
        results = DMA.DMA(betas, samplesheet, group1=GROUP1, group2=GROUP2)
