import os
from collections import defaultdict
import argparse
import re 
import itertools
import numpy as np
from matplotlib import pyplot as plt
from statsmodels.tsa.stattools import acf
from statsmodels.graphics.tsaplots import plot_acf
import pickle
from joblib import Parallel, delayed
import pandas as pd
import seaborn as sns
from datetime import datetime
from multiprocessing import set_start_method, get_context, Manager, Process

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-w", 
        "--bed_file", 
        type=str,
        help="Path to bed file",
        default="None",
        required=False
    )
    parser.add_argument(
        "-wp", 
        "--bed_path", 
        type=str,
        help="Path to folder of multiple bed files",
        default="None",
        required=False
    )
    return parser.parse_args(args)

def parse_bed(filepath: str):
    gene_info = defaultdict()
    current_gene = "current_gene"
    with open(filepath) as fp:
        lines = fp.readlines()
        for line in lines:
            info = re.split('"|\t', line.strip('\n'))
            gene = info[9]
            if gene != current_gene: # If a new gene is being considered
                current_gene = gene
                gene_info[current_gene] = []
            m_val = float(info[15])
            gene_info[gene].append(m_val)
    
    gene_info_by_length = defaultdict()
    for gene, info in gene_info.items():
        exon_count = len(info) 
        if exon_count not in gene_info_by_length:
            gene_info_by_length[exon_count] = [(gene, info)]
        else:
            gene_info_by_length[exon_count].append((gene, info))
    return(gene_info_by_length)
            

def safe_acf(mvals_array: np.array):
    try:
        res = acf(mvals_array)
        return res
    except ValueError:
        pass
    
def compute_acf_parallel(exon_count: str, info: list, acfs_list, genes_list):
    try: 
        gene_ids = [x[0] for x in info]
        genes_list[exon_count] = gene_ids
        acf_res = np.asarray([safe_acf(x[1]) for x in info])
        if acf_res.shape[0] > 1:
            acf_res = acf_res[~np.isnan(acf_res).any(axis=1)]
            acf_res = np.average(acf_res, axis=0)
        else:
            acf_res = acf_res.flatten()
        print("Done for exon length: " + str(exon_count))
            
    except ValueError:
        print("An error occured for exon length: " + str(exon_count))
    acfs_list.append(acf_res)

def get_average_acf(parsed_mvals_dict: dict):
    acfs, genes = [], []
    with Manager() as manager:
        acfs_list = manager.list()
        genes_list = manager.dict()
        processes = []
        for exon_count, info in parsed_mvals_dict.items():
            p = Process(target=compute_acf_parallel, args=(exon_count, info, acfs_list, genes_list))
            p.start()
            processes.append(p)
        for p in processes:
            p.join()
        acfs = list(acfs_list)
        genes = dict(genes_list)
        
    #Sort by exon counts and add padding
    acfs = sorted(acfs, key=len)
    max_exon_count = max(parsed_mvals_dict.keys())
    acfs = [ np.pad(i, (0, max_exon_count - len(i)), 'constant') for i in acfs ]
    acfs = pd.DataFrame(acfs, index=sorted(parsed_mvals_dict.keys())).T
    return (genes, acfs)

def visualize_acf(acfs_by_exon_df: pd.DataFrame):
        acfs_by_exon_long = acfs_by_exon_df.reset_index().melt(id_vars='index', var_name='exon_length', value_name='coef')
        
        sns.set_theme(style="white")
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(14, 10)) 
        
        for i, ax in enumerate(fig.axes):
            custom_exon_count = max(acfs_by_exon_long['exon_length'])/(10**i)
            if i == len(fig.axes)-1:
                custom_exon_count = custom_exon_count*10
                plot_df = acfs_by_exon_long
            else:
                plot_df = acfs_by_exon_long[acfs_by_exon_long['exon_length'] <= custom_exon_count]
            
            plot = sns.lineplot(ax=ax,
                data=plot_df, x='index', y='coef', hue='exon_length',
                palette=sns.color_palette("mako", plot_df['exon_length'].nunique()),
                legend=False, alpha=0.4)
            plot.axhline(0, ls='--', color='red')
            
            
            ax.set(ylabel='Autocorrelation', xlabel='Step',
                    xlim=(0, custom_exon_count), ylim=(-1.2, 1.2),
                    yticks=np.arange(-1,1.2, 0.25))
        
        fig.savefig(file + '_autocor.png', dpi=300)
        plt.clf()
        plt.close(fig)  

def plot_hist_exon_count(freq_list: list):
    fig, ax = plt.subplots(nrows=1, ncols=1) 
    sns.histplot(data=freq_list, ax=ax, palette=sns.color_palette("mako", 1))
    ax.set(xticks=np.arange(1, max(freq_list), 150))
    fig.savefig(file + '_exondist.png', dpi=300)
    plt.clf()
    plt.close(fig)  


    
if __name__ == "__main__":
    set_start_method('spawn')
    args = parse_args()
    bed_path = args.bed_path
    bed_file = args.bed_file
    
    if bed_file == "None" and bed_path == "None":
        print("bed file/path is missing")
        exit()
    elif bed_path != "None":
        print('bed path: ' + bed_path)
        files = [os.path.join(bed_path, f) for f in os.listdir(bed_path) if f.endswith('bed')]
        for file in files:
            print("Started " + file + " at: " + datetime.now().strftime("%H:%M:%S"))  
            # genewise_mval = parse_bed(file)
    
            # acfs_by_exon = get_average_acf(genewise_mval)
            # genes_by_exon = acfs_by_exon[0]
            # acfs_by_exon = acfs_by_exon[1]
            
            # if len(genes_by_exon) > 1:
            #     with open(file + '_genes.pkl', 'wb') as f:
            #         pickle.dump(genes_by_exon, f)
            #     with open(file + '_acfs.pkl', 'wb') as f:
            #         pickle.dump(acfs_by_exon, f)
            
            acfs_by_exon = pickle.load(open(file + '_acfs.pkl', 'rb'))
            # print(np.amax(acfs_by_exon[0:10], axis=1))
            # print(np.mean(acfs_by_exon[0:10], axis=1))
            visualize_acf(acfs_by_exon)
            
            genes_by_exon = pickle.load(open(file + '_genes.pkl', 'rb'))
            freq = [[k]*len(v) for k, v in genes_by_exon.items()]
            freq = list(itertools.chain(*freq))
            plot_hist_exon_count(freq)
            
            
            
    print("Ended: " + datetime.now().strftime("%H:%M:%S"))  