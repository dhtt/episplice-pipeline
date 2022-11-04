import os
from collections import defaultdict
import argparse
import re 
import numpy as np
from statsmodels.graphics.tsaplots import plot_acf

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-w", 
        "--wig_file", 
        type=str,
        help="Path to bigwig/wig file.",
        default="None",
        required=False
    )
    parser.add_argument(
        "-wp", 
        "--wig_path", 
        type=str,
        help="Path to folder of multiple wig files",
        default="None",
        required=False
    )
    return parser.parse_args(args)

def parse_wig(filepath: str):
    chrom = "NA"
    coverage = defaultdict()
    with open(filepath) as fp:
        lines = fp.readlines()
        for line in lines:
            if line.startswith('variableStep'):
                chrom = re.split('=| ', line)[2]
                coverage[chrom] = []
            else: 
                info = re.split('\t', line.strip('\n'))
                coverage[chrom].append(info)
    return(coverage)
            
     

if __name__ == "__main__":
    args = parse_args()
    wig_path = args.wig_path
    wig_file = args.wig_file
    
    if wig_file == "None" and wig_path == "None":
        print("WIG file/path is missing")
        exit()
    elif wig_path != "None":
        print('Wig path: ' + wig_path)
        files = [os.path.join(wig_path, f) for f in os.listdir(wig_path) if f.endswith('wig')]
        for file in files:
            wig_coverage = parse_wig(file)
    
    for chrom, info in enumerate(wig_coverage):
        coverage = np.array(info)
        print(info)
        print(coverage)
        print(coverage.shape)
        # plot_acf()
            