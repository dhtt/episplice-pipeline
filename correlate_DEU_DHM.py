import argparse
from fileinput import filename
import sys
import os
from utils.bash_utils import *
import utils.sam2bed
import pandas as pd
import functools
import operator
from collections import defaultdict
from pathlib import Path
from os import listdir
from os.path import isdir, isfile, join

sys.path.append("/home/dhthutrang/Krebs/episplice-pipeline/utils")

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p1", 
        "--mrna_seq_path", 
        type=str,
        help="Path to mRNA-seq folder."
    )
    parser.add_argument(
        "-p2", 
        "--chip_seq_path", 
        type=str,
        help="Path to ChIP-seq folder."
    )
    parser.add_argument(
        "-gtf1", 
        "--refgen_path", 
        type=str,
        help="Reference genome gtf file"
    )
    parser.add_argument(
        "-gtf2", 
        "--refgen_flank_path", 
        type=str,
        help="Reference genome gtf file contianing only flanks"
    )
    parser.add_argument(
        "-ctl", 
        "--control_id", 
        type=str,
        help="Sample control ID"
    )
    parser.add_argument(
        "-trm", 
        "--treatment_ids", 
        type=str,
        help="Sample treatment IDs"
    )
    return parser.parse_args(args)


def execute_workflow(args=None):
    os.chdir('/home/dhthutrang/Krebs/episplice-pipeline')
    args = parse_args(args)
    print(args)
    
    chip_seq_path = args.chip_seq_path
    mrna_seq_path = args.mrna_seq_path
    control_id = args.control_id
    treatment_ids = args.treatment_ids.split(',')
    refgen_flank_path = args.refgen_flank_path
    refgen_path = args.refgen_path

    cor_analysis_output_path = "/".join(mrna_seq_path.split("/")[:-1]) + "/cor_analysis"
    DHM_output_path = join(chip_seq_path, 'manorm_output_annotated')
    DEU_output_path = '/'.join([mrna_seq_path, 'DEXSeq_output/csv'])

    # fastq_path = join(chip_seq_path, 'fastq_files') #Files in fastq_path must follow naming format celltype_condition_repx_read.fastq.gz
    # chipseq_output = join(chip_seq_path, 'chipseq_output')
    # samplesheet_path = join(fastq_path, 'samplesheet.csv')
    # peak_path = join(chip_seq_path, 'peak_files')
    # alignment_path = join(chip_seq_path, 'bed_files')
    # manorm_output = join(chip_seq_path, 'manorm_output')
    # manorm_output_annotated = join(chip_seq_path, 'manorm_output_annotated')

    # merged_library_path = join(chipseq_output, 'bwa', 'mergedLibrary')
    # narrow_peak_path = join(merged_library_path, 'macs', 'narrowPeak')
    # consensus_peak_path = join(narrow_peak_path, 'consensus')


    print ("========== INITIALIZED DHM WORKFLOW ==========")
    mkdir_p(cor_analysis_output_path)
    histone_types = list()
    for DHM_file in listdir(DHM_output_path):
        histone_types.append(DHM_file.split("_")[0])
    temp = "Rscript get_cor.R --histone_types %s --DEU_file %s --DHM_file %s --output_path %s" % ("_".join(histone_types), DEU_output_path, DHM_output_path, cor_analysis_output_path)
    subprocess.call(temp, shell=True)
   
def test():
    return True

if __name__ == "__main__":
    # test()
    sys.exit(execute_workflow())
