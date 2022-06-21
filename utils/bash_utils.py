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

def merge_bed_file(input_dir: str, output_dir: str, bed_files: list, merged_bed_name: str):
    merged_bed_file_name = join(output_dir, merged_bed_name + '.bed')
    sorted_bed_file_name = join(output_dir, merged_bed_name + '.sorted.bed')
    bed_files = [join(input_dir, s) for s in bed_files]
    print(bed_files)

    subprocess.call(
        "cat %s > %s" %(' '.join(bed_files), merged_bed_file_name),
        shell=True
    )
    subprocess.call(
        "sort -k1,1 -k2,2n %s > %s"
        %(merged_bed_file_name, sorted_bed_file_name),
        shell=True
        )
    subprocess.call(
        "bedtools merge -c 7 -o max -i %s > %s && rm %s"
        %(sorted_bed_file_name, merged_bed_file_name, sorted_bed_file_name),
        shell=True
        )
