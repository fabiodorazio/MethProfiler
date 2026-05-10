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
    parser.add_argument("-q", dest="quality", default=None, help="Lower quality threshold")
    parser.add_argument("--get-mutational-signatures", dest="mut", default=None, help="Run SigProfilerExtractor on vcf input")
    parser.add_argument("--bin", dest="bins", default=20, type=int, help="Number of bins for each chromosome for visualisation purposes")

    args = parser.parse_args()
    beta = args.beta
    samplesheet = args.samplesheet
    name = args.name
    output = args.output_dir
    name = args.name
    test_mode = args.test_mode
    to_keep = args.to_keep
    bin = args.bins

    return (beta, samplesheet, name, output, test_mode, to_keep, bin)


if __name__ == "__main__":
    BETA, SAMPLESHEET, NAME, OUTPUT, TEST_MODE, TO_KEEP, BIN = get_args()

    if NAME is None:
        NAME = utils.get_basename(BETA)

    # checks if output exists or create one
    output_dir = utils.check_output_dir(OUTPUT)
    # creates path for beta
    output_vcf = output_dir + NAME + "_beta.csv"
    # sets location for plot output, checks if exists or creates one
    output_plot_dir = utils.check_output_dir(output_dir + '/Plot_outputs')

    #### MODULES
    betas, samplesheet = DMA.load_meth(BETA,
                SAMPLESHEET,
                test_mode=TEST_MODE,
                samples_to_keep=TO_KEEP)