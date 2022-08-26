import argparse
from email.policy import default
from fileinput import filename
import sys
import os
import subprocess
from utils.bash_utils import *
import utils.sam2bed
import pandas as pd
from collections import defaultdict
from pathlib import Path
from os import listdir
from os.path import isdir, isfile, join

sys.path.append("/home/dhthutrang/Krebs/episplice-pipeline/utils")
import fastq_dir_to_samplesheet #noLint


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p", 
        "--chip_seq_path", 
        type=str,
        help="Path to mRNA-seq folder."
    )
    parser.add_argument(
        "-r1", 
        "--read1_extension", 
        type=str,
        help="File extension for read 1."
    )
    parser.add_argument(
        "-r2", 
        "--read2_extension", 
        type=str,
        help="File extension for read 2."
    )
    parser.add_argument(
        "-G", 
        "--genome", 
        type=str,
        help="Reference genome version"
    )
    parser.add_argument(
        "-ctl", 
        "--control_id", 
        type=str,
        help="Sample control ID"
    )
    parser.add_argument(
        "-trm", 
        "--treatment_id", 
        type=str,
        help="Sample treatment IDs"
    )
    parser.add_argument(
        "-ab", 
        "--antibody", 
        type=str,
        help="ChIP-seq antibody"
    )
    parser.add_argument(
        "-ri",
        "--replicate_index",
        type=int,
        default=1,
        help="After splitting FastQ file name by --sanitise_name_delimiter all elements before this index (1-based) will be joined to create final sample name.",
    )
    parser.add_argument(
        "-gi",
        "--group_name_index",
        type=str,
        default="0",
        help="After splitting FastQ file name by --sanitise_name_delimiter all elements before this index (1-based) will be joined to create final sample name.",
    )
    parser.add_argument(
        "-nco",
        "--chipseq_options",
        type=str,
        default="0",
        help="After splitting FastQ file name by --sanitise_name_delimiter all elements before this index (1-based) will be joined to create final sample name.",
    )
    return parser.parse_args(args)


def execute_workflow(args=None):
    os.chdir('/home/dhthutrang/Krebs/episplice-pipeline')
    args = parse_args(args)
    print(args)
    
    chip_seq_path = args.chip_seq_path
    read1_extension = args.read1_extension
    read2_extension = args.read2_extension
    genome = args.genome
    control_id = args.control_id
    treatment_id = args.treatment_id.split(',')
    antibody = str(args.antibody)
    replicate_index = str(args.replicate_index)
    group_name_index = str(args.group_name_index)
    chipseq_options = args.chipseq_options
    
    # #TEST VALUES
    # chip_seq_path = '/home/dhthutrang/Krebs/Messier/ChIP_seq'
    # narrow_peak_path = '/home/dhthutrang/Krebs/Reddy/data/ChIP_seq/chipseq_output/bwa/mergedLibrary/macs/broadPeak/consensus'
    # read1_extension = '.fastq.gz'
    # read2_extension = '.fastq.gz'
    # genome = 'GRCh38'
    # control_id = 'MCF7_ctrl'
    # treatment_id = ['MCF7_e2_', 'MCF7_gc10_']
    # antibody = 2
    # replicate_index = 3
    # group_name_index = '0,1,2'
    
    fastq_path = join(chip_seq_path, 'fastq_files') #Files in fastq_path must follow naming format celltype_condition_repx_read.fastq.gz
    chipseq_output = join(chip_seq_path, 'chipseq_output')
    samplesheet_path = join(fastq_path, 'samplesheet.csv')
    peak_path = join(chip_seq_path, 'peak_files')
    bed_path = join(chip_seq_path, 'bed_files')

    merged_library_path = join(chipseq_output, 'bwa', 'mergedLibrary')
    narrow_peak_path = join(merged_library_path, 'macs', 'narrowPeak')
    consensus_peak_path = join(narrow_peak_path, 'consensus')


    # print ("========== Initialized DHM workflow ==========")
    # # Run nf-core/chip_seq
    # print("Create output folder")
    # mkdir_p(chipseq_output)
    
    # print("Make sample sheet")
    # fastq_dir_to_samplesheet.main(args = [fastq_path, samplesheet_path, 'chipseq',
    #     '--read1_extension', read1_extension, 
    #     '--read2_extension', read2_extension,
    #     '--control_id', control_id,
    #     '--antibody', antibody,
    #     '--replicate_index', replicate_index,
    #     '--group_name_index', group_name_index
    #     ])
        
    # print("Run nf-core/chipseq")
    # subprocess.call(
    #     "nextflow run nf-core/chipseq --input %s --outdir %s --genome %s -profile docker --save_reference true"
    #     %(samplesheet_path, chipseq_output, genome) +
    #     " " + chipseq_options, 
    #     shell=True)
    

    # Prepare MAnorm
    print("Create MAnorm input folder")
    mkdir_p(peak_path)
    mkdir_p(bed_path)
    histone_types = set()
    all_consensus_file_id = set()
    for sorted_sam_file in listdir(merged_library_path):
        if sorted_sam_file.endswith('sorted.bam'):
            narrow_peak_file_name = sorted_sam_file.split('.')[0]
            histone_types.add(narrow_peak_file_name.split("_")[int(antibody)])
            all_consensus_file_id.add('_'.join(narrow_peak_file_name.split('_')[:-1]) )


    print("PREPARE PEAK BED FILES")
    # For each control sample, peak file is generated by merging narrow peak files 
    print("1. Control sample from merged consensus peaks")
    for consensus_file_id in all_consensus_file_id: 
        narrow_peak_files = [file for file in listdir(narrow_peak_path) if file.startswith(consensus_file_id) & file.endswith('peaks.narrowPeak') & (control_id in file)]
        narrow_peak_files = [file for file in narrow_peak_files if len(file)>0]
        if (len(narrow_peak_files) != 0):
            merge_bed_file(narrow_peak_path, peak_path, narrow_peak_files, consensus_file_id)


    # For each treatment sample, peak file is generated by merging MACS consensus peak files 
    print("2. Treatment sample from MACS consensus peaks")
    for histone_type in histone_types:
        if histone_type not in listdir(consensus_peak_path): 
            file = [f for f in listdir(narrow_peak_path) if (f.endswith(".narrowPeak")) & (histone_type in f) & (control_id not in f)][0]
            mkdir_p(join(consensus_peak_path, histone_type))
            duplicate_file(join(narrow_peak_path, file), join(consensus_peak_path, histone_type, histone_type+'.consensus_peaks.boolean.txt'))
    
        file = join(consensus_peak_path, histone_type, histone_type + '.consensus_peaks.boolean.txt')
        # Read consensus peak file and get sample names
        df = pd.read_csv(file, sep='\t', header=0)
        sample_name = [col for col in df.columns if col.endswith('.bool')]
        sample_name = sample_name[0] if len(sample_name) > 1 else df.iloc[0, 3]
        sample_name = sample_name.split('_')
        print(sample_name)
        if len(sample_name) > 4: 
            file_name = '_'.join(sample_name[:-3]) + '.bed'
        else:
            file_name = '_'.join(sample_name[:-1]) + '.bed'
        (df
        .iloc[:, :3]
        .to_csv(join(peak_path, file_name), sep='\t', header=False, index=False)
        )

    # print("PREPARE READ ALIGNMENT BED FILES")
    # print("Converting reads SAM to BED")
    # # Convert read alignments from BAM to SAM then SAM to BED
    # for file in listdir(merged_library_path):
    #     if file.endswith('sorted.bam'):
    #         file_name = file.split('.')[0]
    #         file_name_sam = file_name + ".sam"
    #         file_name_bed = file_name + ".bed"

    #         subprocess.call(
    #             "./utils/bam2sam.sh -i %s -o %s"
    #             %(merged_library_path, bed_path),
    #             shell=True
    #         )
    #         subprocess.call(
    #             "python utils/sam2bed.py -i %s -o %s"
    #             %(join(bed_path, file_name_sam), join(bed_path, file_name_bed)),
    #             shell=True
    #         )
            
    # # Merge BEDs of different replicate in one BED
    # for histone_type in histone_types:
    #     sample_ids = treatment_id + [control_id]
    #     for sample_id in sample_ids:
    #         file_id = '_'.join([sample_id, histone_type])    
    #         alignment_files = [file for file in listdir(bed_path) if file.startswith(file_id) & file.endswith('.bed')]
    #         merge_bed_file(bed_path, bed_path, alignment_files, file_id)

    print("========== Finished nf-core/chipseq ==========")

def test():
    return True

if __name__ == "__main__":
    # test()
    sys.exit(execute_workflow())
