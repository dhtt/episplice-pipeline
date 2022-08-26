import subprocess
from pathlib import Path
from os import listdir
from os.path import isdir, isfile, join

def mkdir_p(path_to_create):
    Path(path_to_create).mkdir(parents=True, exist_ok=True)


def duplicate_file(current_path, new_path):
    subprocess.call(
        "cp %s %s" %(current_path, new_path),
        shell=True
        )

def sort_bed_file(bed_file):
    subprocess.call(
        "sort -k1,1 -k2,2n %s > %s"
        %(bed_file, bed_file + ".sorted"),
        shell=True
        )
    subprocess.call(
        "mv %s %s"
        %(bed_file + ".sorted", bed_file),
        shell=True
        )

def rm_file(file):
    subprocess.call(
        "rm %s" %(file),
        shell=True
        )

def merge_bed_file(input_dir: str, output_dir: str, bed_files: list, merged_bed_name: str):
    # Define input and outname
    merged_bed_file_name = join(output_dir, merged_bed_name + '.bed')
    bed_files = [join(input_dir, s) for s in bed_files]

    # Join bed files by cat
    subprocess.call(
        "cat %s > %s" %(' '.join(bed_files), merged_bed_file_name),
        shell=True
    )
    
    # Sort bed files by bedtools
    sort_bed_file(merged_bed_file_name)
    subprocess.call(
        "bedtools merge -c 5 -o max -i %s > %s"
        %(merged_bed_file_name, merged_bed_file_name + ".merged"),
        shell=True
        )
    subprocess.call(
        "mv %s %s"
        %(merged_bed_file_name + ".merged", merged_bed_file_name),
        shell=True
        )
