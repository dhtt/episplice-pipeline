import argparse
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


    print ("========== Initialized DHM workflow ==========")
    # Run nf-core/chip_seq
    print("Create output folder")
    mkdir_p(chipseq_output)
    
    print("Make sample sheet")
    fastq_dir_to_samplesheet.main(args = [fastq_path, samplesheet_path, 'chipseq',
        '--read1_extension', read1_extension, 
        '--read2_extension', read2_extension,
        '--control_id', control_id,
        '--antibody', antibody,
        '--replicate_index', replicate_index,
        '--group_name_index', group_name_index
        ])
        
    print("Run nf-core/chipseq")
    subprocess.call(
        "nextflow run nf-core/chipseq --input %s --outdir %s --genome %s -profile docker --save_reference true"
        %(samplesheet_path, chipseq_output, genome) +
        " " + chipseq_options, 
        shell=True)
    

    # Prepare MAnorm
    print("Create MAnorm input folder")
    mkdir_p(peak_path)
    mkdir_p(bed_path)


    print("Call peaks for control sample")
    control_consensus_file_id = set()
    for sorted_sam_file in listdir(merged_library_path):
        if sorted_sam_file.startswith(control_id) & sorted_sam_file.endswith('sorted.bam'):
            narrow_peak_file_name = sorted_sam_file.split('.')[0]

            # If sample are single-end, remove only replication index
            # consensus_file_name is often celltype_condition_antibody
            if read2_extension == read1_extension:
                consensus_file_name = '_'.join(narrow_peak_file_name.split('_')[:-1]) 
            # else if sample are paired-end, remove replication index and read index
            else: 
                consensus_file_name = '_'.join(narrow_peak_file_name.split('_')[:-2]) 
            control_consensus_file_id.add(consensus_file_name)

            BAM_type = "BAM" if read1_extension==read2_extension else "BAMPE"
            sorted_sam_path = join(merged_library_path, sorted_sam_file)
            subprocess.call(
                "macs2 callpeak -t %s -c %s -f %s -g 2.7e9 -n %s --outdir %s --keep-dup all"
                %(sorted_sam_path, sorted_sam_path, BAM_type, narrow_peak_file_name, narrow_peak_path), 
                shell=True
                )


    print("Merge consensus peaks for control sample")
    for consensus_file_id in control_consensus_file_id: 
        narrow_peak_files = [file for file in listdir(narrow_peak_path) if file.startswith(consensus_file_id) & file.endswith('peaks.narrowPeak')]
        merge_bed_file(narrow_peak_path, peak_path, narrow_peak_files, consensus_file_id)

    print("Filtering MACS peaks")
    # Process consensus peak files from consensus_peak_path to peak_path for each antibody
    for histone_type in listdir(consensus_peak_path):
        file = join(consensus_peak_path, histone_type, histone_type + '.consensus_peaks.boolean.txt')
        # Read consensus peak file and get sample names
        df = pd.read_csv(file, sep='\t', header=0)
        samples = [col for col in df.columns if col.endswith('.bool')]

        for s in samples:
            # If sample are single-end, remove only replication index
            if read2_extension == read1_extension:
                file_name = '_'.join(s.split('_')[:-1]) + '.bed' 
            # else if sample are paired-end, remove replication index and read index
            else: 
                file_name = '_'.join(s.split('_')[:-2]) + '.bed' 

            # Write chr, start, end to BED file    
            df[df[s] == True].iloc[:, :3].to_csv(
                join(peak_path, file_name),
                sep='\t', header=False, index=False
            )

    print("Converting reads SAM to BED")
    # Convert read alignments from BAM to SAM then SAM to BED
    for file in listdir(merged_library_path):
        if file.endswith('sorted.bam'):
            file_name = file.split('.')[0]
            file_name_sam = file_name + ".sam"
            
            # If sample are single-end, remove only replication index
            if read2_extension == read1_extension:
                file_name_bed = '_'.join(file_name.split('_')[:-1]) + '.bed' 
            # else if sample are paired-end, remove replication index and read index
            else: 
                file_name_bed = '_'.join(file_name.split('_')[:-2]) + '.bed'

            subprocess.call(
                "./utils/bam2sam.sh -i %s -o %s"
                %(merged_library_path, bed_path),
                shell=True
            )
            subprocess.call(
                "python utils/sam2bed.py -i %s -o %s"
                %(join(bed_path, file_name_sam), join(bed_path, file_name_bed)),
                shell=True
            )
           
    print("========== Finished nf-core/chipseq ==========")

def test():
    return True

if __name__ == "__main__":
    # test()
    sys.exit(execute_workflow())
